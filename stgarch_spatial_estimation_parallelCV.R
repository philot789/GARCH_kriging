# STGARCH(1,1): local QMLE and the two feasible predictors in Appendix C.
# Only complete, finite input data are supported. Preprocessing is external.
# Coordinates must be on the same two-dimensional scale used in the paper.
# The xi covariance below uses the working relation C_xi = 2 C_eta^2;
# that relation is an additional fourth-moment restriction, not distribution-free.

stgarch_bandwidth <- function(dx, dy) {
  sqrt(mean(c(stats::var(dx), stats::var(dy)))) * length(dx)^(-1 / 6)
}

# backward compatibility
default_bandwidth <- function(dx, dy) {
  coords <- cbind(dx, dy)
  sigma <- sqrt(mean(apply(coords, 2, var)))
  
  sigma * length(dx)^(-1 / 6)
}

# backward compatibility
a2d_sepgauss <- function(dx, dy, b) {
  log_weights <- -0.5 * (
    (dx / b)^2 +
      (dy / b)^2
  )
  
  weights <- exp(log_weights - max(log_weights))
  
  weights / sum(weights)
}

# backward compatibility
# Define uniform (box) kernels
a2d_box <- function(dx, dy, b) {
  # rectangular window: |dx| ≤ b and |dy| ≤ b
  w <- (abs(dx) <= b) & (abs(dy) <= b)
  w / sum(w)
}
# backward compatibility
# (optional circular window)
a2d_circle <- function(dx, dy, b) {
  d2 <- dx^2 + dy^2
  w <- (d2 <= b^2)
  w / sum(w)
}
# backward compatibility
# gaussian 
a2d_gaussian <- function(dx, dy, b) {
  
  d2 <- dx^2 + dy^2
  w  <- exp( - d2 / (2 * b^2) )
  w / sum(w)
}


# GARCH simulation given eta and params
simulation <- function(eta, omega, alpha, beta) {
  T <- nrow(eta);
  n <- ncol(eta)
  
  Z <- sigma2 <- matrix(0, nrow = T, ncol = n)
  
  sigma2[1, ] <- omega
  
  for (t in 2:T) {
    for (u in 1:n) {
      sigma2[t,u] <- omega[u] + alpha[u]*Z[t-1,u]^2 + beta[u]*sigma2[t-1,u]
      Z[t,u] <- sqrt(sigma2[t,u]) * eta[t,u]
    }
  }
  return(list(Z = Z, sigma2 = sigma2))
}




stgarch_weights <- function(dx, dy, b) {
  # A separable biweight kernel: symmetric, Lipschitz, and compactly supported,
  # as required by the theoretical kernel assumptions. Normalizing the weights
  # does not change the local QMLE minimizer.
  w <- pmax(1 - (dx / b)^2, 0)^2 * pmax(1 - (dy / b)^2, 0)^2
  if (sum(w) > 0) w / sum(w) else w
}

.stgarch_data <- function(Z, coords) {
  Z <- as.matrix(Z)
  coords <- as.matrix(coords)
  if (!is.numeric(Z) || !is.numeric(coords) || nrow(Z) < 3L ||
      ncol(Z) < 2L || nrow(coords) != ncol(Z) || ncol(coords) != 2L ||
      any(!is.finite(Z)) || any(!is.finite(coords)) ||
      any(!is.finite(Z^2))) {
    stop("Supply finite numeric Z (time x sites) and coords (sites x 2), ",
         "with at least three times and two sites.")
  }
  if (anyDuplicated(as.data.frame(coords)))
    stop("coords must identify distinct sites; duplicate coordinates need an explicit nugget interpretation.")
  list(Z = Z, coords = coords)
}

.stgarch_controls <- function(theta_init, lower, upper, stationarity_margin) {
  if (!is.numeric(theta_init) || !is.numeric(lower) || !is.numeric(upper) ||
      length(theta_init) != 3L || length(lower) != 3L || length(upper) != 3L ||
      any(!is.finite(c(theta_init, lower, upper))) ||
      lower[1] <= 0 || any(lower[2:3] < 0) || any(upper <= lower) ||
      length(stationarity_margin) != 1L || !is.finite(stationarity_margin) ||
      stationarity_margin <= 0 || stationarity_margin >= 1 ||
      any(theta_init <= lower) || any(theta_init >= upper) ||
      sum(theta_init[2:3]) >= 1 - stationarity_margin) {
    stop("Use finite bounds, lower[1] > 0, nonnegative alpha/beta bounds, ",
         "and an interior theta_init with alpha + beta < 1 - stationarity_margin.")
  }
}

.stgarch_parallel <- function(parallel, n_workers) {
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel))
    stop("parallel must be TRUE or FALSE.")
  if (!parallel) return(function() invisible(NULL))
  if (!requireNamespace("future", quietly = TRUE) ||
      !requireNamespace("future.apply", quietly = TRUE))
    stop("Install the future and future.apply packages for parallel fitting.")
  if (length(n_workers) != 1L || !is.finite(n_workers) ||
      n_workers < 2L || n_workers != as.integer(n_workers))
    stop("Use at least two integer workers when parallel = TRUE.")
  old_plan <- future::plan("list")
  # Separate R sessions work on Windows, macOS, and Linux, including RStudio.
  future::plan(future::multisession, workers = as.integer(n_workers))
  function() future::plan(old_plan)
}

.stgarch_map <- function(X, FUN, parallel) {
  if (parallel) {
    future.apply::future_lapply(X, FUN, future.seed = TRUE)
  } else {
    lapply(X, FUN)
  }
}

.stgarch_bandwidth <- function(coords, bandwidth, bandwidth_function) {
  if (is.null(bandwidth)) bandwidth <- bandwidth_function(coords[, 1], coords[, 2])
  if (length(bandwidth) != 1L || !is.finite(bandwidth) || bandwidth <= 0)
    stop("The bandwidth must be a finite positive scalar.")
  bandwidth
}

.stgarch_qml <- function(theta, Z2, w) {
  # The same unconditional initialization is used for fitting and inversion.
  # Row 1 initializes the recursion; rows 2:T enter the likelihood.
  omega <- theta[1]; alpha <- theta[2]; beta <- theta[3]
  gap <- 1 - alpha - beta
  if (omega <= 0 || alpha < 0 || beta < 0 || gap <= 0)
    return(list(value = Inf, gradient = rep(NA_real_, 3)))
  h <- rep(omega / gap, ncol(Z2))
  dh <- matrix(rep(c(1 / gap, omega / gap^2, omega / gap^2),
                   each = ncol(Z2)), ncol = 3L)
  value <- 0
  gradient <- numeric(3)
  for (tt in 2:nrow(Z2)) {
    dh <- cbind(1, Z2[tt - 1L, ], h) + beta * dh
    h <- omega + alpha * Z2[tt - 1L, ] + beta * h
    value <- value + sum(w * (log(h) + Z2[tt, ] / h)) / 2
    gradient <- gradient + colSums(dh * (w * (1 - Z2[tt, ] / h) / (2 * h)))
  }
  list(value = value / (nrow(Z2) - 1L),
       gradient = gradient / (nrow(Z2) - 1L))
}

.stgarch_local <- function(Z2, coords, target, theta_init, lower, upper,
                           bandwidth, weight_function, stationarity_margin) {
  w <- weight_function(coords[, 1] - target[1], coords[, 2] - target[2], bandwidth)
  if (length(w) != ncol(Z2) || any(!is.finite(w)) || any(w < 0) || sum(w) <= 0)
    stop("Invalid spatial weights or no sites inside the kernel support.")
  w <- w / sum(w)
  # Scaling changes the optimization coordinates, not the model or objective.
  scale <- c(max(sum(colMeans(Z2) * w), theta_init[1]), 1, 1)
  ui <- rbind(diag(3), -diag(3), c(0, -1, -1))
  ci <- c(lower, -upper, -(1 - stationarity_margin))
  ui <- sweep(ui, 2L, scale, `*`)
  last_x <- NULL; last_result <- NULL
  evaluate <- function(x) {
    if (!identical(x, last_x)) {
      last_result <<- .stgarch_qml(x * scale, Z2, w)
      last_x <<- x
    }
    last_result
  }
  fit <- stats::constrOptim(
    theta = theta_init / scale,
    f = function(x) evaluate(x)$value,
    grad = function(x) evaluate(x)$gradient * scale,
    ui = ui, ci = ci, method = "BFGS",
    control = list(maxit = 1000L, reltol = 1e-9),
    outer.iterations = 100L, outer.eps = 1e-7
  )
  theta <- fit$par * scale
  if (fit$convergence != 0L || !is.finite(fit$value) ||
      any(!is.finite(theta)) || any(theta < lower) || any(theta > upper) ||
      sum(theta[2:3]) > 1 - stationarity_margin)
    stop("Local QMLE failed at (", paste(signif(target, 6), collapse = ", "),
         "). Optimizer code: ", fit$convergence, ". No prediction was substituted.")
  stats::setNames(theta, c("omega", "alpha", "beta"))
}

stgarch_invert <- function(Z, parameters) {
  # parameters must correspond, in order, to the observed columns of Z.
  Z <- as.matrix(Z); parameters <- as.matrix(parameters)
  if (!is.numeric(Z) || !is.numeric(parameters) || nrow(Z) < 2L ||
      nrow(parameters) != ncol(Z) || ncol(parameters) != 3L ||
      any(!is.finite(Z^2)) || any(!is.finite(parameters)) ||
      any(parameters[, 1] <= 0) || any(parameters[, 2:3] < 0) ||
      any(rowSums(parameters[, 2:3, drop = FALSE]) >= 1))
    stop("Invalid data or parameters for inversion.")
  h <- matrix(NA_real_, nrow(Z), ncol(Z), dimnames = dimnames(Z))
  h[1, ] <- parameters[, 1] / (1 - parameters[, 2] - parameters[, 3])
  eta <- matrix(NA_real_, nrow(Z), ncol(Z))
  for (tt in 2:nrow(Z)) {
    h[tt, ] <- parameters[, 1] + parameters[, 2] * Z[tt - 1L, ]^2 +
      parameters[, 3] * h[tt - 1L, ]
    eta[tt, ] <- Z[tt, ] / sqrt(h[tt, ])
  }
  if (any(!is.finite(h)) || any(!is.finite(eta[-1L, , drop = FALSE])))
    stop("Nonfinite variance or standardized residuals.")
  list(sigma2 = h, eta = eta, xi = eta^2 - 1)
}

.stgarch_covariance <- function(eta, D) {
  # Fit the exponential-plus-nugget model in (27) by Gaussian quasi-likelihood.
  # The free total variance is a nuisance scale; the fitted covariance is then
  # normalized to C_eta(0) = 1. Residual values are NOT rescaled or re-centered.
  positive_D <- D[D > 0]
  if (!length(positive_D)) stop("Spatial covariance requires distinct coordinates.")
  scale_D <- stats::median(positive_D)
  D_scaled <- D / scale_D
  S <- crossprod(eta) / nrow(eta)
  marginal <- mean(diag(S))
  if (!is.finite(marginal) || marginal <= 0) stop("Degenerate standardized residuals.")
  objective <- function(par) {
    R <- par[2] * exp(-D_scaled / par[1]) + diag(par[3], nrow(D))
    U <- tryCatch(chol(R), error = function(e) NULL)
    if (is.null(U)) return(.Machine$double.xmax / 100)
    inv_R <- chol2inv(U)
    sum(log(diag(U))) + sum(inv_R * S) / 2
  }
  fit <- stats::optim(
    par = c(1, 0.8 * marginal, 0.2 * marginal), fn = objective,
    method = "L-BFGS-B", lower = c(min(positive_D) / scale_D / 100, 1e-10, 1e-10),
    upper = c(2 * max(D_scaled), Inf, Inf),
    control = list(maxit = 1000L, factr = 1e7)
  )
  if (fit$convergence != 0L || !is.finite(fit$value) || any(!is.finite(fit$par)))
    stop("Spatial covariance fit failed; optimizer code: ", fit$convergence)
  range <- fit$par[1] * scale_D
  fraction <- fit$par[2] / sum(fit$par[2:3])
  R_eta <- fraction * exp(-D / range) + diag(1 - fraction, nrow(D))
  if (fit$par[1] <= min(positive_D) / scale_D / 100 * (1 + 1e-6) ||
      fit$par[1] >= 2 * max(D_scaled) * (1 - 1e-6))
    warning("Spatial range is at a search boundary; inspect covariance adequacy.")
  list(range = range, spatial_variance = fraction, nugget_variance = 1 - fraction,
       raw_total_variance = sum(fit$par[2:3]), R_eta = R_eta)
}

.stgarch_predict <- function(theta, eta, r_eta, R_chol, predictor, epsilon) {
  r <- if (predictor == "xi_blp") 2 * r_eta^2 else r_eta
  gamma <- as.numeric(backsolve(R_chol, forwardsolve(t(R_chol), r)))
  marginal <- if (predictor == "xi_blp") 2 else 1
  v <- marginal - sum(r * gamma)
  if (!is.finite(v) || v < -1e-8 * marginal)
    stop("Invalid kriging variance; check the covariance matrix and target covariance.")
  v <- max(v, 0) # Only remove negative roundoff; not a positivity adjustment.
  projection <- as.numeric(
    (if (predictor == "xi_blp") eta[-1L, , drop = FALSE]^2 - 1 else
      eta[-1L, , drop = FALSE]) %*% gamma
  )
  h <- z2 <- numeric(nrow(eta))
  spatial_mspe <- rep(NA_real_, nrow(eta))
  h[1] <- z2[1] <- theta[1] / (1 - theta[2] - theta[3])
  for (tt in 2:nrow(eta)) {
    h[tt] <- theta[1] + theta[2] * z2[tt - 1L] + theta[3] * h[tt - 1L]
    multiplier <- if (predictor == "xi_blp") 1 + projection[tt - 1L] else
      projection[tt - 1L]^2 + v
    raw <- h[tt] * multiplier
    if (!is.finite(h[tt]) || h[tt] <= 0 || !is.finite(raw))
      stop("Nonfinite prediction; no floor or fallback can repair this fit.")
    # Equation (35) truncates the PRODUCT, not the standardized multiplier.
    z2[tt] <- if (predictor == "xi_blp") max(raw, epsilon) else raw
    # (36) is for the unadjusted spatial error component, not total uncertainty.
    # Gaussian expression follows Var(eta_0^2 | eta) = 2 v^2 + 4 mu^2 v.
    spatial_mspe[tt] <- h[tt]^2 * if (predictor == "xi_blp") v else
      (2 * v^2 + 4 * projection[tt - 1L]^2 * v)
  }
  if (any(z2 < 0) || any(!is.finite(spatial_mspe[-1L])))
    stop("Invalid prediction or spatial error variance.")
  list(z2 = z2, h = h, spatial_mspe = spatial_mspe)
}

fit_stgarch_spatial_cv <- function(
    Z, coords, theta_init, lower, upper,
    predictor = c("xi_blp", "eta_gaussian"),
    folds = NULL, n_folds = 5L, seed = 5515L,
    weight_function = stgarch_weights, bandwidth_function = stgarch_bandwidth,
    bandwidth = NULL, stationarity_margin = 1e-4, epsilon = 1e-8,
    parallel = FALSE, n_workers = 2L, progress = interactive()) {
  predictor <- match.arg(predictor)
  input <- .stgarch_data(Z, coords); Z <- input$Z; coords <- input$coords
  .stgarch_controls(theta_init, lower, upper, stationarity_margin)
  if (length(epsilon) != 1L || !is.finite(epsilon) || epsilon <= 0)
    stop("epsilon must be finite and positive.")
  n <- ncol(Z)
  if (is.null(folds)) {
    if (length(n_folds) != 1L || !is.finite(n_folds) || n_folds < 2L ||
        n_folds > n || n_folds != as.integer(n_folds)) stop("Invalid n_folds.")
    # Preserve the caller's random-number state when generating folds.
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    folds <- sample(rep(seq_len(n_folds), length.out = n))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else
      rm(".Random.seed", envir = .GlobalEnv)
  }
  if (length(folds) != n || anyNA(folds) || length(unique(folds)) < 2L)
    stop("folds must assign every site to one of at least two folds.")
  bandwidth <- .stgarch_bandwidth(coords, bandwidth, bandwidth_function)
  restore <- .stgarch_parallel(parallel, n_workers); on.exit(restore(), add = TRUE)
  Z_pred_S2 <- sigma2_pred_S2 <- spatial_mspe <- matrix(NA_real_, nrow(Z), n,
                                                        dimnames = dimnames(Z))
  parameters <- matrix(NA_real_, n, 3L,
                       dimnames = list(colnames(Z), c("omega", "alpha", "beta")))
  covariance_fits <- list()
  for (fold in unique(folds)) {
    if (progress) message("Fitting fold ", fold)
    test <- which(folds == fold); train <- which(folds != fold)
    if (length(train) < 2L) stop("Every fold needs at least two training sites.")
    Z_train <- Z[, train, drop = FALSE]
    coords_train <- coords[train, , drop = FALSE]
    local_fit <- function(target) .stgarch_local(
      Z_train^2, coords_train, target, theta_init, lower, upper,
      bandwidth, weight_function, stationarity_margin
    )
    training_parameters <- do.call(rbind, .stgarch_map(
      seq_along(train), function(j) local_fit(coords_train[j, ]), parallel
    ))
    eta <- stgarch_invert(Z_train, training_parameters)$eta
    D <- as.matrix(stats::dist(coords_train))
    cov_fit <- .stgarch_covariance(eta[-1L, , drop = FALSE], D)
    R <- if (predictor == "xi_blp") 2 * cov_fit$R_eta^2 else cov_fit$R_eta
    R_chol <- chol(R)
    predictions <- .stgarch_map(test, function(i) {
      theta <- local_fit(coords[i, ])
      distances <- sqrt(rowSums(sweep(coords_train, 2L, coords[i, ], `-`)^2))
      # No nugget in cross-covariance to a distinct, unobserved site.
      r_eta <- cov_fit$spatial_variance * exp(-distances / cov_fit$range)
      list(theta = theta, prediction = .stgarch_predict(
        theta, eta, r_eta, R_chol, predictor, epsilon
      ))
    }, parallel)
    for (j in seq_along(test)) {
      i <- test[j]; out <- predictions[[j]]
      parameters[i, ] <- out$theta
      Z_pred_S2[, i] <- out$prediction$z2
      sigma2_pred_S2[, i] <- out$prediction$h
      spatial_mspe[, i] <- out$prediction$spatial_mspe
    }
    cov_fit$R_eta <- NULL
    covariance_fits[[as.character(fold)]] <- cov_fit
  }
  # The empirical section defines the volatility score using sqrt(h), not
  # sqrt(predicted Z^2). These are different prediction targets. Absolute Z is
  # a noisy proxy, so this is not direct error against the latent true sigma.
  error <- sqrt(sigma2_pred_S2[-1L, , drop = FALSE]) - abs(Z[-1L, , drop = FALSE])
  if (any(!is.finite(error))) stop("Incomplete CV evaluation; refusing a partial score.")
  if (!is.finite(mean(error^2))) stop("Overflow in the CV score.")
  list(Z_pred_S2 = Z_pred_S2, sigma2_pred_S2 = sigma2_pred_S2,
       estimated_params = parameters, estimated_params_all = parameters,
       spatial_mspe = spatial_mspe, folds = folds, predictor = predictor,
       RMSE_vol = sqrt(mean(error^2)), MAE_vol = mean(abs(error)),
       evaluation_coverage = 1, bandwidth = bandwidth, covariance_fits = covariance_fits)
}

fit_stgarch_local_parameters <- function(
    Z, coords, theta_init, lower, upper, target_coords = coords,
    weight_function = stgarch_weights, bandwidth_function = stgarch_bandwidth,
    bandwidth = NULL, stationarity_margin = 1e-4,
    parallel = FALSE, n_workers = 2L, progress = interactive()) {
  input <- .stgarch_data(Z, coords); Z <- input$Z; coords <- input$coords
  .stgarch_controls(theta_init, lower, upper, stationarity_margin)
  target_coords <- as.matrix(target_coords)
  if (!is.numeric(target_coords) || ncol(target_coords) != 2L ||
      nrow(target_coords) < 1L || any(!is.finite(target_coords)))
    stop("target_coords must be a finite numeric matrix with two columns.")
  bandwidth <- .stgarch_bandwidth(coords, bandwidth, bandwidth_function)
  restore <- .stgarch_parallel(parallel, n_workers); on.exit(restore(), add = TRUE)
  if (progress) message("Fitting ", nrow(target_coords), " local parameter vectors")
  parameters <- do.call(rbind, .stgarch_map(seq_len(nrow(target_coords)), function(i) {
    .stgarch_local(Z^2, coords, target_coords[i, ], theta_init, lower, upper,
                   bandwidth, weight_function, stationarity_margin)
  }, parallel))
  list(estimated_params = parameters, estimated_params_all = parameters,
       target_coords = target_coords, bandwidth = bandwidth)
}
