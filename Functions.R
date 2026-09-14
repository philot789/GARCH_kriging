# Rule-of-thumb bandwidth for 2D data (Scott's rule)
default_bandwidth <- function(dx, dy) {
  coords <- cbind(dx, dy)
  sigma <- sqrt(mean(apply(coords, 2, var)))  # pooled std dev
  n <- length(dx)
  sigma * n^(-1/6)
}

#separable‐Gaussian weight
a2d_sepgauss <- function(dx, dy, b) {
  # dx = x_u - x0, dy = y_u - y0, b = window‐width
  kx <- exp(- (dx/b)^2 / 2)
  ky <- exp(- (dy/b)^2 / 2)
  w  <- kx * ky
  w / sum(w)
}

# Define uniform (box) kernels
a2d_box <- function(dx, dy, b) {
  # rectangular window: |dx| ≤ b and |dy| ≤ b
  w <- (abs(dx) <= b) & (abs(dy) <= b)
  w / sum(w)
}
# (optional circular window)
a2d_circle <- function(dx, dy, b) {
  d2 <- dx^2 + dy^2
  w <- (d2 <= b^2)
  w / sum(w)
}

#gaussian 
a2d_gaussian <- function(dx, dy, b) {
  
  d2 <- dx^2 + dy^2
  w  <- exp( - d2 / (2 * b^2) )
  w / sum(w)
}






# Likelihood computation
compute_likelihood <- function(theta, Z, weights) {
  omega <- theta[1];
  alpha <- theta[2]; 
  beta <- theta[3]
  T <- nrow(Z); 
  m <- ncol(Z)
  
  sigma2 <- matrix(0, nrow = T, ncol = m)
  sigma2[1, ] <- omega
  
  p <- 1
  
  for (t in (p + 1):T) {
    sigma2[t, ] <- omega + alpha * Z[t-1, ]^2 + beta * sigma2[t-1, ]
  }
  
  ll_terms <- 0.5 * (log(sigma2[(p + 1):T, ]) + (Z[(p + 1):T, ]^2) / sigma2[(p + 1):T, ])
  
  weighted_ll <- sum(t(t(ll_terms) * as.numeric(weights))) / ((T - p - 1) * m)
  
  return(weighted_ll)
}

# Inversion to get eta_hat and zeta_hat
inversion_eta_hat <- function(Z, omega, alpha, beta) {
  T <- nrow(Z);
  n <- ncol(Z)
  
  sigma2 <- eta_hat <- zeta_hat <- matrix(0, nrow = T, ncol = n)
  sigma2[1, ] <- omega
  
  for (t in 2:T) {
    for (u in 1:n) {
      sigma2[t,u] <- omega[u] + alpha[u]*Z[t-1,u]^2 + beta[u]*sigma2[t-1,u]
      eta_hat[t,u] <- Z[t,u] / sqrt(sigma2[t,u])
      zeta_hat[t,u] <- sigma2[t,u] * (eta_hat[t,u]^2 - 1)
    }
  }
  
  eta_hat[1, ] <- NA
  zeta_hat[1, ] <- NA
  
  return(list(eta = eta_hat, zeta = zeta_hat))
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

# Full negative log-likelihood for spatial range
fullNegLogLik <- function(theta_range, D, eta) {
  tau2 <- 1e-3;
  n <- ncol(eta); 
  T <- nrow(eta)
  C <- exp(-D / theta_range) + tau2 * diag(n)
  L <- chol(C);
  Cinv <- chol2inv(L)
  logdet <- 2 * sum(log(diag(L)))
  Q <- sum(apply(eta,1,function(y) t(y) %*% Cinv %*% y))
  0.5 * T * logdet + 0.5 * Q
}

# Full negative log-likelihood for spatial range and variance
fullNegLogLik_sig_theta <- function(theta, D, eta) {
  theta_range <- theta[1]
  sigma2 <- theta[2]
  tau2 <- 1e-3
  n <- ncol(eta)
  T <- nrow(eta)
  
  C <- sigma2 * exp(-D / theta_range) + tau2 * diag(n)
  L <- chol(C)
  Cinv <- chol2inv(L)
  logdet <- 2 * sum(log(diag(L)))
  Q <- sum(apply(eta, 1, function(y) t(y) %*% Cinv %*% y))
  
  0.5 * T * logdet + 0.5 * Q
}

fullNegLogLik_sig_theta_tau <- function(theta, D, eta) {
  # Extract parameters
  theta_range <- theta[1]
  sigma2 <- theta[2]
  tau2 <- theta[3]
  
  # Sanity checks for parameter bounds
  if (any(theta <= 0)) return(1e10)  # Enforce positivity
  
  n <- ncol(eta)
  T <- nrow(eta)
  
  # Construct covariance matrix
  C <- sigma2 * exp(-D / theta_range) + tau2 * diag(n)
  
  # Cholesky decomposition with error handling
  L <- tryCatch(chol(C), error = function(e) return(NULL))
  if (is.null(L)) return(1e10)
  
  # Quadratic form using backsolve (more stable than inversion)
  Q <- sum(apply(eta, 1, function(y) sum(backsolve(L, y, transpose = TRUE)^2)))
  
  # Log determinant
  logdet <- 2 * sum(log(diag(L)))
  
  # Negative log-likelihood
  0.5 * T * logdet + 0.5 * Q
}


compute_likelihood_log <- function(theta, Z, weights) {
  omega <- theta[1]
  alpha <- theta[2]
  beta  <- theta[3]
  
  T <- nrow(Z)
  m <- ncol(Z)
  
  log_sigma2 <- matrix(0, nrow = T, ncol = m)
  log_sigma2[1, ] <- omega
  
  for (t in 2:T) {
    log_sigma2[t, ] <- omega +
      alpha * log(Z[t - 1, ]^2 + 1e-8) +
      beta  * log_sigma2[t - 1, ]
  }
  
  ll_terms <- 0.5 * (log_sigma2[2:T, ] + Z[2:T, ]^2 / exp(log_sigma2[2:T, ]))
  weighted_ll <- sum(t(t(ll_terms) * weights)) / ((T - 1) * m)
  
  return(weighted_ll)  # to be minimised
}


inversion_eta_hat_log <- function(Z, omega, alpha, beta) {
  T <- nrow(Z)
  n <- ncol(Z)
  
  log_sigma2 <- matrix(0, nrow = T, ncol = n)
  eta_hat <- matrix(0, nrow = T, ncol = n)
  zeta_hat <- matrix(0, nrow = T, ncol = n)
  
  log_sigma2[1, ] <- omega  # consistent with simulation
  
  for (t in 2:T) {
    for (u in 1:n) {
      log_sigma2[t, u] <- omega[u] +
        alpha[u] * log(Z[t - 1, u]^2 + 1e-8) +
        beta[u]  * log_sigma2[t - 1, u]
      eta_hat[t, u] <- Z[t, u] / sqrt(exp(log_sigma2[t, u]))
      zeta_hat[t, u] <- exp(log_sigma2[t, u]) * (eta_hat[t, u]^2 - 1);
    }
  }
  eta_hat[1, ] <- NA
  zeta_hat[1, ] <- NA
  
  return(list(eta = eta_hat, zeta = zeta_hat))
}

simulation_log <- function(eta, omega, alpha, beta) {
  
  T <- nrow(eta)
  n <- ncol(eta)
  
  log_sigma2 <- matrix(0, nrow = T, ncol = n)
  Z <- matrix(0, nrow = T, ncol = n)
  
  log_sigma2[1, ] <- omega  # initial log-variance
  
  for (t in 2:T) {
    for (u in 1:n) {
      log_sigma2[t, u] <- omega[u] +
        alpha[u] * log(Z[t - 1, u]^2 + 1e-8) +  # note that we assume that P(Z_true = 0) = 0, thus we add a little bit to ensure that Z_true is non-zero
        beta[u]  * log_sigma2[t - 1, u]
      Z[t, u] <- sqrt(exp(log_sigma2[t, u])) * eta[t, u]
    }
  }
  
  return(list(Z = Z, log_sigma2 = log_sigma2))
}




# Tensor-product construction: row-wise Kronecker products
kronecker_basis <- function(Bx, By) {
  n <- nrow(Bx);
  p1 <- ncol(Bx); 
  p2 <- ncol(By)
  B  <- matrix(0, nrow=n, ncol=p1*p2)
  for(i in 1:n) {
    B[i,] <- as.vector(outer(Bx[i,], By[i,])) 
  }
  return(B)
}