#loading required packages
library(stats)
library(parallel)
library(pbapply)
library(spGARCH)
library(Matrix)
library(Rsolnp)
library(rlist)


weight_matrix<-function(covariates){
  covariates=as.matrix(covariates)
  n=nrow(covariates)
  weight_matrix=matrix(0,ncol = n,nrow=n)
  for(i in 1:n){
    for(j in 1:n){
      weight_matrix[i,j] = 1/sqrt(sum((covariates[i,] - covariates[j,])^2))
      
      if(i==j){
        weight_matrix[i,j]=0
      }
    }
    
  }
  #standardizing
  
  #weight_matrix=(1/max(weight_matrix[weight_matrix!=Inf]))*weight_matrix
  # weight_matrix[weight_matrix==Inf]=1
  weight_matrix=weight_matrix*array(rep(1/colSums(weight_matrix),n),dim = c(n,n))
  
  return(weight_matrix)
}


weight_matrix<-function(covariates){
  covariates=as.matrix(covariates)
  n=nrow(covariates)
  weight_matrix=matrix(0,ncol = n,nrow=n)
  for(i in 1:n){
    for(j in 1:n){
      weight_matrix[i,j] = 1/sqrt(sum((covariates[i,] - covariates[j,])^2))
      
      if(i==j){
        weight_matrix[i,j]=0
      }
    }
    
  }
  #weight_matrix[weight_matrix==Inf]=max(weight_matrix[weight_matrix!=Inf])
  #standardizing
  
  weight_matrix[weight_matrix==Inf]=1e12
  
  
  weight_matrix=diag(1/rowSums(weight_matrix))%*%weight_matrix
  #weight_matrix=(1/max(weight_matrix[weight_matrix!=Inf]))*weight_matrix
  #weight_matrix[weight_matrix==Inf]=1
  rowSums(weight_matrix)
  return(weight_matrix)
}


choose_functions <- function(model){
  
  if (model == "spGARCH")
  {
    
    f_inv <- function(x){
      return(x)
    }
    tau_eps <- function(eps, W_1, W_2, alpha, theta, zeta, b){  # function tau in paper
      eps_2 <- eps^2
      n     <- length(eps)
      return(diag(eps_2) %*% solve(diag(n) - Matrix(W_1) %*% diag(eps_2) - Matrix(W_2)) %*% alpha)
    }
    tau_y <- function(y){
      return(y^2)
    }
    g <- function(y, alpha, rhoW_1, lambdaW_2, theta, zeta, b, tau_y, f_inv){ # mapping from y to eps
      n   <- length(y)
      X   <- solve(diag(n) - Matrix(lambdaW_2)) %*% (alpha + Matrix(rhoW_1) %*% tau_y(y))
      eps <- y / sqrt(as.vector(f_inv(X)))
      return(list(eps = eps, h = f_inv(X)))
    }
    d_h_d_eps <- function(eps, h, alpha, rhoW_1, lambdaW_2){
      n <- length(eps)
      aux <- solve(diag(n) - rhoW_1 %*% diag(as.vector(eps^2)) - lambdaW_2)
      result <- array(, dim = c(n,n))
      for(j in 1:n){
        aux_rhoW_1 <- rhoW_1
        aux_rhoW_1[,-j] <- array(0, dim = c(n, n-1))
        eps_j <- eps[j]
        result[, j] <- 2 * eps_j * aux %*% aux_rhoW_1 %*% aux %*% alpha
      }
      return(result)
    }

  } 
  else if (model == "spEGARCH")
  {
    
    f_inv <- function(x){
      return(exp(x))
    }
    tau_eps <- function(eps, W_1, W_2, alpha, theta, zeta, b){  # function tau in paper
      return(log(abs(eps)^b))
    }
    tau_y <- function(y){
      return(NULL)
    }
    g <- function(y, alpha, rhoW_1, lambdaW_2, theta, zeta, b, tau_y, f_inv){ # mapping from y to eps
      # is only valid if g(eps) = (ln(|eps_1|^b), ..., ln(|eps_n|^b))
      n     <- length(y)
      # b     <- 2 
      # X     <- solve(diag(n) + Matrix(rhoW_1) - Matrix(lambdaW_2)) %*% (alpha + Matrix(rhoW_1) %*% log(abs(y)^b))
      X     <- solve(diag(n) + 0.5 * b * Matrix(rhoW_1) - Matrix(lambdaW_2)) %*% (alpha + b * Matrix(rhoW_1) %*% log(abs(y)))
      eps   <- y / sqrt(as.vector(f_inv(X)))
      return(list(eps = eps, h = f_inv(X)))
    }
    d_h_d_eps <- function(eps, h, alpha, theta, zeta, b, rhoW_1, lambdaW_2){
      tau_eps_prime <- function(eps){
        return(b / eps) # Achtung evtl. zu ändern, hard coded!
      }
      n     <- length(eps)
      aux   <- solve(diag(n) - Matrix(lambdaW_2)) %*% Matrix(rhoW_1)
      eps_j <- t(array(rep(eps, n), dim = c(n, n)))
      h_i   <- array(rep(h, n), dim = c(n, n))
      return(aux * tau_eps_prime(eps_j) * h_i)
    }
    
  }
  else if (model == "spEGARCH2")
  {
    
    f_inv <- function(x){
      return(exp(x))
    }
    tau_eps <- function(eps, W_1, W_2, alpha, theta, zeta, b){  # function tau in paper
      return(theta * eps)
    }
    tau_y <- function(y){
      return(NULL)
    }
    g <- function(y, alpha, rhoW_1, lambdaW_2, theta, zeta, b, tau_y, f_inv){ # mapping from y to eps
      n     <- length(y)
      # X     <- solve(diag(n) + Matrix(rhoW_1) - Matrix(lambdaW_2)) %*% (alpha + Matrix(rhoW_1) %*% log(abs(y)^b))
      # mapping from y to eps (or rather X = log(eps))
      function_y_eps <- function(x, y, alpha, rhoW_1, lambdaW_2, theta, zeta){
        n <- length(x)
        eps <- x
        return(    (diag(n) - lambdaW_2) %*% log(abs(y)^2) - alpha     -     ( rhoW_1 %*% (theta * eps) + (diag(n) - lambdaW_2) %*% log(abs(eps)^2) )  )
      }
      # function_y_eps(x = eps, y = Y, alpha = alpha, rhoW_1 = rhoW_1, lambdaW_2 = lambdaW_2)
      out <- nleqslv(x = y, fn = function_y_eps, alpha = alpha, rhoW_1 = rhoW_1, lambdaW_2 = lambdaW_2, theta = theta, zeta = zeta, y = y, control = list(ftol = 1e-10))
      # plot(eps, out$x)
      # plot((out$x), ylim = c(-4,4))
      # points((eps), col = "red")
      # X     <- solve(diag(n) + 0.5 * b * Matrix(rhoW_1) - Matrix(lambdaW_2)) %*% (alpha + b * Matrix(rhoW_1) %*% log(abs(y)))
      # eps   <- y / sqrt(as.vector(f_inv(X)))
      eps <- out$x
      X <- log(abs(y)^2) - log(abs(eps)^2)
      return(list(eps = eps, h =  f_inv(X)))
    }
    d_h_d_eps <- function(eps, h, alpha, theta, zeta, b, rhoW_1, lambdaW_2){
      tau_eps_prime <- function(eps, theta){
        return(theta) # Achtung evtl. zu ändern, hard coded!
      }
      n     <- length(eps)
      aux   <- solve(diag(n) - Matrix(lambdaW_2)) %*% Matrix(rhoW_1)
      eps_j <- t(array(rep(eps, n), dim = c(n, n)))
      h_i   <- array(rep(h, n), dim = c(n, n))
      return(aux * tau_eps_prime(eps_j, theta) * h_i)
    }
    
  }
  else if (model == "spHARCH")
  {
    
    f_inv <- function(x){
      return(exp(x))
    }
    tau_eps <- function(eps, W_1, W_2, alpha, theta, zeta, b){  # function tau in paper
      eps_2 <- log(eps^2)
      n     <- length(eps)
      return( solve(diag(n) - Matrix(W_1) - Matrix(W_2)) %*% alpha + (diag(n) + solve(diag(n) - Matrix(W_1) - Matrix(W_2)) %*% Matrix(W_1)) %*% eps_2  )
    } 
    tau_y <- function(y){
      log(y^2)
    }
    g <- function(y, alpha, rhoW_1, lambdaW_2, theta, zeta, b, tau_y, f_inv){ # mapping from y to eps
      n   <- length(y)
      X   <- solve(diag(n) - Matrix(lambdaW_2)) %*% (alpha + rhoW_1 %*% tau_y(y))
      eps <- y / sqrt(as.vector(f_inv(X)))
      return(list(eps = eps, h = f_inv(X)))
    }
    d_h_d_eps <- function(eps, h, alpha, theta, zeta, b, rhoW_1, lambdaW_2){
      n     <- length(eps)
      aux   <- solve(diag(n) - Matrix(rhoW_1) - Matrix(lambdaW_2)) %*% Matrix(rhoW_1)
      eps_j <- t(array(rep(eps, n), dim = c(n, n)))
      h_i   <- array(rep(h, n), dim = c(n, n))
      return(2 * aux / eps_j * h_i)
    }
  }
  
  return(list(f_inv = f_inv, tau_eps = tau_eps, tau_y = tau_y, g = g, d_h_d_eps = d_h_d_eps))
}



model="spHARCH"

functions <- choose_functions(model)
f_inv     <- functions$f_inv
g         <- functions$g
tau_eps   <- functions$tau_eps
tau_y     <- functions$tau_y
d_h_d_eps <- functions$d_h_d_eps

if(is.element(model, c("spGARCH", "spHARCH"))){
  pars     <- c(0.5, 0.5, 1)
  LB       <- c(0, 0, 0.00001)
  UB       <- c(Inf, Inf, Inf)
} else if(is.element(model, c("spEGARCH"))) {
  pars     <- c(0.5, 0.5, 1)
  LB       <- c(0, 0, 0.00001)
  UB       <- c(Inf, Inf, Inf)
} else if(is.element(model, c("spEGARCH2"))) {
  pars     <- c(0.5, 0.5, 1, 0.5)
  LB       <- c(0, 0, 0.00001, 0.00001)
  UB       <- c(Inf, Inf, Inf, Inf)
}

Least_Squares <- function(pars_K, param){
  
  rho    <- pars_K[1]
  lambda <- pars_K[2]
  alpha  <- pars_K[3]
  theta  <- 1
  zeta   <- 1
  b      <- 1
  
  y         <- param$y
  f_inv     <- param$f_inv
  tau_y     <- param$tau_y
  W_1       <- param$W_1
  W_2       <- param$W_2
  g         <- param$g
  
  n            <- length(y)
  rhoW_1       <- rho * W_1
  lambdaW_2    <- lambda * W_2
  
  result_g     <- g(y, alpha, rhoW_1, lambdaW_2, theta, zeta, b, tau_y, f_inv)
  h            <- as.vector(result_g[[2]])
  
  z <- log(y^2) - (digamma(1) - log(2))
  h_tilde <- log(h)
  
  squares <- (z - h_tilde)^2    
  
  return(sum(squares))
  
}

LL <- function(pars, param){
  
  model  <- param$model
  
  if(is.element(model, c("spGARCH", "spHARCH"))){
    rho    <- pars[1]
    lambda <- pars[2]
    alpha  <- pars[3]
    theta  <- 1
    zeta   <- 1
    b      <- 1
  } else if(is.element(model, c("spEGARCH"))) {
    rho    <- pars[1]
    lambda <- pars[2]
    alpha  <- pars[3]
    theta  <- NULL
    zeta   <- NULL
    b      <- pars[4]
  } else if(is.element(model, c("spEGARCH2"))) {
    rho    <- pars[1]
    lambda <- pars[2]
    alpha  <- pars[3]
    theta  <- pars[4]
    zeta   <- 0 # pars[5]
    b      <- NULL
  }
  
  
  y         <- param$y
  f_inv     <- param$f_inv
  tau_y     <- param$tau_y
  W_1       <- param$W_1
  W_2       <- param$W_2
  g         <- param$g
  d_h_d_eps <- param$d_h_d_eps
  
  n            <- length(y)
  rhoW_1       <- rho * W_1
  lambdaW_2    <- lambda * W_2
  alpha        <- alpha * rep(1, n)
  
  result_g     <- g(y, alpha, rhoW_1, lambdaW_2, theta, zeta, b, tau_y, f_inv)
  eps          <- as.vector(result_g[[1]])
  h            <- as.vector(result_g[[2]])
  
  J <- 0.5 * eps / sqrt(h) * d_h_d_eps(eps, h, alpha, theta, zeta, b, rhoW_1, lambdaW_2) + diag(sqrt(h))
  
  # print(pars)
  
  log_det_J <- determinant(J, logarithm = TRUE)$modulus  
  
  return(as.numeric((-1) * (sum(log(dnorm(eps))) -  log_det_J)))
  
}

#set number of replications
m = 1000



#set number of objects/firms to be located
n = 50

sim_wrapper <- function(i,n, omega, rho, lambda, alpha){
  #compute true locations
  S = matrix(runif(n*2,0,1),nrow = n, ncol = 2)
  WM_true = weight_matrix(S)
  

  #simulate series
  
  eps   <- rnorm(n)
  
  tau   <- tau_eps(eps, rho * WM_true, lambda * WM_true, rep(alpha, n), theta, zeta, b)
  Inv   <- solve(diag(n) - lambda * WM_true)
  X     <- Inv %*% (alpha + rho * WM_true %*% tau)
  
  sim_series     <- diag(as.vector(f_inv(X)))^(0.5) %*% eps
  
  # some tests
  # result_g     <- g(sim_series, alpha, rho * WM_true, lambda * WM_true, theta, zeta, b, tau_y, f_inv)
  # h            <- as.vector(result_g[[2]])
  # round(log(h), 5) == round(as.vector(X), 5)
  # round(as.vector(result_g[[1]]), 5) == round(eps, 5)
  
  
  # to plot simulated process:
  # summary(sim_series)
  # colfunc <- colorRampPalette(c(rgb(1,1,0.7), "orange", rgb(0.6, 0, 0)))
  # cols <- colfunc(30)
  # q8 <- classIntervals(sim_series, n = 30, style = "fisher")
  # q8Cols <- findColours(q8,cols)
  # 
  # plot(S, col = q8Cols, pch = 20, lwd = 10)
  
  #sim_series = spGARCH::sim.spGARCH(n = n, rho = rho, lambda = lambda, alpha = alpha, W1 = Matrix(WM_true), W2 = Matrix(WM_true), type = "log-spGARCH")
  #sim_series = spGARCH::sim.spARCH(n = n, rho = rho, alpha = alpha, W = WM_true, type = "spARCH")
  
  X_star = (1-omega)*S + omega * matrix(rnorm(n*2),nrow = n, ncol = 2)
  
  #getting R^2
  df_lm = data.frame(S,X_star)
  r1 = lm(X1~X1.1, df_lm)
  r2 = lm(X2~X2.1, df_lm)
  r1=summary(r1)$r.squared 
  r2=summary(r2)$r.squared 
  WM_estimate = weight_matrix(X_star)
  
  #estimating parameter
  param    <- list(y = sim_series, f_inv = functions$f_inv, tau_y = functions$tau_y, d_h_d_eps = functions$d_h_d_eps, g = functions$g, W_1 = WM_estimate, W_2 = WM_estimate, model = model)
  out      <- solnp(c(rho, lambda, alpha), fun = LL, LB = c(0, 0, 0.00001), UB = c(1, 1, 100), param = param, control = list(trace = FALSE))
  
  # out$pars / sqrt(diag(solve(out$hessian)))
  # round(out$pars, 3)
  
  res = c(out$pars,r1,r2)
  return(res)
}
#set number of cores
nc=10

cl <- makeCluster(nc)
setDefaultCluster(cl)
#setting seed for reproducibility

set.seed(4572)
sim_results_setup1 = pblapply(FUN=sim_wrapper, X= 1:m, n=50, omega=0, rho=0.7, lambda=0.2, alpha=1)
sim_results_setup2 = pblapply(FUN=sim_wrapper, X= 1:m, n=100, omega=0, rho=0.7, lambda=0.2, alpha=1)
sim_results_setup3 = pblapply(FUN=sim_wrapper, X= 1:m, n=150, omega=0, rho=0.7, lambda=0.2, alpha=1)

sim_results_setup1=list.rbind(sim_results_setup1)
sim_results_setup2=list.rbind(sim_results_setup2)
sim_results_setup3=list.rbind(sim_results_setup3)

set.seed(4572)
sim_results_setup1_50_omega_01 = pblapply(FUN=sim_wrapper, X= 1:m, n=50, omega=0.1, rho=0.7, lambda=0.2, alpha=1)
sim_results_setup1_100_omega_01 = pblapply(FUN=sim_wrapper, X= 1:m, n=100, omega=0.1, rho=0.7, lambda=0.2, alpha=1)
sim_results_setup1_150_omega_01 = pblapply(FUN=sim_wrapper, X= 1:m, n=150, omega=0.1, rho=0.7, lambda=0.2, alpha=1)

sim_results_setup1_50_omega_01=list.rbind(sim_results_setup1_50_omega_01)
sim_results_setup1_100_omega_01=list.rbind(sim_results_setup1_100_omega_01)
sim_results_setup1_150_omega_01=list.rbind(sim_results_setup1_150_omega_01)

set.seed(4572)
sim_results_setup1_50_omega_02 = pblapply(FUN=sim_wrapper, X= 1:m, n=50, omega=0.2, rho=0.7, lambda=0.2, alpha=1)
sim_results_setup1_100_omega_02 = pblapply(FUN=sim_wrapper, X= 1:m, n=100, omega=0.2, rho=0.7, lambda=0.2, alpha=1)
sim_results_setup1_150_omega_02 = pblapply(FUN=sim_wrapper, X= 1:m, n=150, omega=0.2, rho=0.7, lambda=0.2, alpha=1)

set.seed(4572)
sim_results_setup1_50_omega_05 = pblapply(FUN=sim_wrapper, X= 1:m, n=50, omega=0.5, rho=0.7, lambda=0.2, alpha=1)
sim_results_setup1_100_omega_05 = pblapply(FUN=sim_wrapper, X= 1:m, n=100, omega=0.5, rho=0.7, lambda=0.2, alpha=1)
sim_results_setup1_150_omega_05 = pblapply(FUN=sim_wrapper, X= 1:m, n=150, omega=0.5, rho=0.7, lambda=0.2, alpha=1)

sim_results_setup1_50_omega_05=list.rbind(sim_results_setup1_50_omega_05)
sim_results_setup1_100_omega_05=list.rbind(sim_results_setup1_100_omega_05)
sim_results_setup1_150_omega_05=list.rbind(sim_results_setup1_150_omega_05)
save(sim_results_setup1_50_omega_05,file="sim_results_setup1_50_omega_05.RData")
save(sim_results_setup1_100_omega_05,file="sim_results_setup1_100_omega_05.RData")
save(sim_results_setup1_150_omega_05,file="sim_results_setup1_150_omega_05.RData")

save(sim_results_setup1,file="sim_results_setup1.RData")
save(sim_results_setup2,file="sim_results_setup2.RData")
save(sim_results_setup3,file="sim_results_setup3.RData")

stopCluster(cl)

set.seed(4572)
sim_results_setup2_50_omega_05 = pblapply(FUN=sim_wrapper, X= 1:m, n=50, omega=0.5, rho=0.3, lambda=0.4, alpha=0.5)
sim_results_setup2_100_omega_05 = pblapply(FUN=sim_wrapper, X= 1:m, n=100, omega=0.5, rho=0.3, lambda=0.4, alpha=0.5)
sim_results_setup2_150_omega_05 = pblapply(FUN=sim_wrapper, X= 1:m, n=150, omega=0.5, rho=0.3, lambda=0.4, alpha=0.5)
sim_results_setup2_50_omega_05=list.rbind(sim_results_setup2_50_omega_05)
sim_results_setup2_100_omega_05=list.rbind(sim_results_setup2_100_omega_05)
sim_results_setup2_150_omega_05=list.rbind(sim_results_setup2_150_omega_05)
save(sim_results_setup2_50_omega_05,file="sim_results_setup2_50_omega_05.RData")
save(sim_results_setup2_100_omega_05,file="sim_results_setup2_100_omega_05.RData")
save(sim_results_setup2_150_omega_05,file="sim_results_setup2_150_omega_05.RData")


#mean
round(mean(sim_results_setup1[,1]),3)
round(mean(sim_results_setup1[,2]),3)
round(mean(sim_results_setup1[,3]),3)

#sd
round(sd(sim_results_setup1[,1]),3)
round(sd(sim_results_setup1[,2]),3)
round(sd(sim_results_setup1[,3]),3)


#RMSE
round(sqrt(mean((sim_results_setup1[,1]-0.7)^2)),3)
round(sqrt(mean((sim_results_setup1[,2]-0.2)^2)),3)
round(sqrt(mean((sim_results_setup1[,3]-1)^2)),3)

#mean
round(mean(sim_results_setup2[,1]),3)
round(mean(sim_results_setup2[,2]),3)
round(mean(sim_results_setup2[,3]),3)

#sd
round(sd(sim_results_setup2[,1]),3)
round(sd(sim_results_setup2[,2]),3)
round(sd(sim_results_setup2[,3]),3)


#RMSE
round(sqrt(mean((sim_results_setup2[,1]-0.7)^2)),3)
round(sqrt(mean((sim_results_setup2[,2]-0.2)^2)),3)
round(sqrt(mean((sim_results_setup2[,3]-1)^2)),3)

#mean
round(mean(sim_results_setup3[,1]),3)
round(mean(sim_results_setup3[,2]),3)
round(mean(sim_results_setup3[,3]),3)

#sd
round(sd(sim_results_setup3[,1]),3)
round(sd(sim_results_setup3[,2]),3)
round(sd(sim_results_setup3[,3]),3)


#RMSE
round(sqrt(mean((sim_results_setup3[,1]-0.7)^2)),3)
round(sqrt(mean((sim_results_setup3[,2]-0.2)^2)),3)
round(sqrt(mean((sim_results_setup3[,3]-1)^2)),3)



#0mega = 0.1
#mean
round(mean(sim_results_setup1_50_omega_01[,1]),3)
round(mean(sim_results_setup1_50_omega_01[,2]),3)
round(mean(sim_results_setup1_50_omega_01[,3]),3)

#sd
round(sd(sim_results_setup1[,1]),3)
round(sd(sim_results_setup1[,2]),3)
round(sd(sim_results_setup1[,3]),3)


#RMSE
round(sqrt(mean((sim_results_setup1_50_omega_01[,1]-0.7)^2)),3)
round(sqrt(mean((sim_results_setup1_50_omega_01[,2]-0.2)^2)),3)
round(sqrt(mean((sim_results_setup1_50_omega_01[,3]-1)^2)),3)

#mean
round(mean(sim_results_setup1_100_omega_01[,1]),3)
round(mean(sim_results_setup1_100_omega_01[,2]),3)
round(mean(sim_results_setup1_100_omega_01[,3]),3)

#sd
round(sd(sim_results_setup1_100_omega_01[,1]),3)
round(sd(sim_results_setup1_100_omega_01[,2]),3)
round(sd(sim_results_setup1_100_omega_01[,3]),3)


#RMSE
round(sqrt(mean((sim_results_setup1_100_omega_01[,1]-0.7)^2)),3)
round(sqrt(mean((sim_results_setup1_100_omega_01[,2]-0.2)^2)),3)
round(sqrt(mean((sim_results_setup1_100_omega_01[,3]-1)^2)),3)

#mean
round(mean(sim_results_setup1_150_omega_01[,1]),3)
round(mean(sim_results_setup1_150_omega_01[,2]),3)
round(mean(sim_results_setup1_150_omega_01[,3]),3)

#sd
round(sd(sim_results_setup1_150_omega_01[,1]),3)
round(sd(sim_results_setup1_150_omega_01[,2]),3)
round(sd(sim_results_setup1_150_omega_01[,3]),3)


#RMSE
round(sqrt(mean((sim_results_setup1_150_omega_01[,1]-0.7)^2)),3)
round(sqrt(mean((sim_results_setup1_150_omega_01[,2]-0.2)^2)),3)
round(sqrt(mean((sim_results_setup1_150_omega_01[,3]-1)^2)),3)


#omega=0.5
#mean
round(mean(sim_results_setup1_50_omega_05[,1]),3)
round(mean(sim_results_setup1_50_omega_05[,2]),3)
round(mean(sim_results_setup1_50_omega_05[,3]),3)

#sd
round(sd(sim_results_setup1_50_omega_05[,1]),3)
round(sd(sim_results_setup1_50_omega_05[,2]),3)
round(sd(sim_results_setup1_50_omega_05[,3]),3)


#RMSE
round(sqrt(mean((sim_results_setup1_50_omega_05[,1]-0.7)^2)),3)
round(sqrt(mean((sim_results_setup1_50_omega_05[,2]-0.2)^2)),3)
round(sqrt(mean((sim_results_setup1_50_omega_05[,3]-1)^2)),3)

#mean
round(mean(sim_results_setup1_100_omega_05[,1]),3)
round(mean(sim_results_setup1_100_omega_05[,2]),3)
round(mean(sim_results_setup1_100_omega_05[,3]),3)

#sd
round(sd(sim_results_setup1_100_omega_05[,1]),3)
round(sd(sim_results_setup1_100_omega_05[,2]),3)
round(sd(sim_results_setup1_100_omega_05[,3]),3)


#RMSE
round(sqrt(mean((sim_results_setup1_100_omega_05[,1]-0.7)^2)),3)
round(sqrt(mean((sim_results_setup1_100_omega_05[,2]-0.2)^2)),3)
round(sqrt(mean((sim_results_setup1_100_omega_05[,3]-1)^2)),3)


#mean
round(mean(sim_results_setup1_150_omega_05[,1]),3)
round(mean(sim_results_setup1_150_omega_05[,2]),3)
round(mean(sim_results_setup1_150_omega_05[,3]),3)

#sd
round(sd(sim_results_setup1_150_omega_05[,1]),3)
round(sd(sim_results_setup1_150_omega_05[,2]),3)
round(sd(sim_results_setup1_150_omega_05[,3]),3)


#RMSE
round(sqrt(mean((sim_results_setup1_150_omega_05[,1]-0.7)^2)),3)
round(sqrt(mean((sim_results_setup1_150_omega_05[,2]-0.2)^2)),3)
round(sqrt(mean((sim_results_setup1_150_omega_05[,3]-1)^2)),3)


#mean
round(mean(sim_results_setup2_50_omega_05[,1]),3)
round(mean(sim_results_setup2_50_omega_05[,2]),3)
round(mean(sim_results_setup2_50_omega_05[,3]),3)

#sd
round(sd(sim_results_setup2_50_omega_05[,1]),3)
round(sd(sim_results_setup2_50_omega_05[,2]),3)
round(sd(sim_results_setup2_50_omega_05[,3]),3)


#RMSE
round(sqrt(mean((sim_results_setup2_50_omega_05[,1]-0.3)^2)),3)
round(sqrt(mean((sim_results_setup2_50_omega_05[,2]-0.4)^2)),3)
round(sqrt(mean((sim_results_setup2_50_omega_05[,3]-0.5)^2)),3)

#mean
round(mean(sim_results_setup2_100_omega_05[,1]),3)
round(mean(sim_results_setup2_100_omega_05[,2]),3)
round(mean(sim_results_setup2_100_omega_05[,3]),3)

#sd
round(sd(sim_results_setup2_100_omega_05[,1]),3)
round(sd(sim_results_setup2_100_omega_05[,2]),3)
round(sd(sim_results_setup2_100_omega_05[,3]),3)


#RMSE
round(sqrt(mean((sim_results_setup2_100_omega_05[,1]-0.3)^2)),3)
round(sqrt(mean((sim_results_setup2_100_omega_05[,2]-0.4)^2)),3)
round(sqrt(mean((sim_results_setup2_100_omega_05[,3]-0.5)^2)),3)


#mean
round(mean(sim_results_setup2_150_omega_05[,1]),3)
round(mean(sim_results_setup2_150_omega_05[,2]),3)
round(mean(sim_results_setup2_150_omega_05[,3]),3)

#sd
round(sd(sim_results_setup2_150_omega_05[,1]),3)
round(sd(sim_results_setup2_150_omega_05[,2]),3)
round(sd(sim_results_setup2_150_omega_05[,3]),3)


#RMSE
round(sqrt(mean((sim_results_setup2_150_omega_05[,1]-0.3)^2)),3)
round(sqrt(mean((sim_results_setup2_150_omega_05[,2]-0.4)^2)),3)
round(sqrt(mean((sim_results_setup2_150_omega_05[,3]-0.5)^2)),3)
