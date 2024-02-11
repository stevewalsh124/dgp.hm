
# Function contents -----------------------------------------------------------
# External:
#   bohman
#   matrix.Moore.Penrose
#   matix.Moore.Penrose2
#   matrix.sqrt

# Optimize Matern -------------------------------------------------------------
#' @export
#' 
opt_matern <- function(dx, y, sdd, init = c(0.1, -5, 0.1, 0.1), 
                       lower = c(-5, -35, -5, 0.2),
                       upper = c(4, 3, 4, 100)) { 
  out <- optim(init, nl_matern, method = "L-BFGS-B", lower = lower,
               upper = upper, dx = dx, y = y, sdd = sdd)
  return(list(kappa_hat = out$par[4], 
              phi_hat = exp(out$par[1]),
              theta_hat = exp(out$par[1]) * out$par[4], 
              g_hat = exp(out$par[2]),
              tau2_hat = exp(out$par[3])))
}

# Negative log likelihood Matern ----------------------------------------------

nl_matern <- function(par, dx, y, sdd) {
  # par: c(theta, g, tau2, kappa)
  # dx: squared distances of inputs
  # y: vector of response
  # sdd: scaled precisions
  
  theta <- exp(par[1])
  g <- exp(par[2])
  tau2 <- exp(par[3])
  kappa <- par[4]
  n <- nrow(y)
  d <- ncol(y)
  
  K <- tau2 * (geoR::matern(sqrt(dx), phi = theta, kappa = kappa) + diag(g, n))
  Ki <- solve(K)
  ldetK <- determinant(K, logarithm = TRUE)$modulus
  
  ll <- 0
  for (i in 1:d) {
    yit <- y[, i] / sdd
    ll <- ll - (1/2)*(t(yit) %*% Ki %*% yit) - (1/2)*ldetK 
  }
  return(-ll)
}

# scrP ------------------------------------------------------------------------
#' @export

scrP <- function(x, k) {
  y <- log10(x * (k^1.5) / (2 * pi^2))
  y[is.infinite(y)] <- 0
  return(y)
}

# bohman ----------------------------------------------------------------------
#' @export

bohman <- function(t, tau = 0.25) {
  boh <- (1 - t / tau) * cos(pi * t / tau) + sin(pi * t / tau) / pi
  boh <- ifelse(t >= tau, 0, boh)
  return(boh)
}

# matrix.Moore.Penrose --------------------------------------------------------
# Computes Moore-Penrose generalized inverse of symmetric matrix using 
# spectral decomposition
# By M.A.R. Ferreira
#' @export

matrix.Moore.Penrose <- function(H) {
  H.eigen <- eigen(H)
  inverse.values <- rep(0, nrow(H))
  inverse.values[abs(H.eigen$values) > 10^(-10)] <- 1 / H.eigen$values[abs(H.eigen$values) > 10^(-10)]
  H.MP <- H.eigen$vectors %*% diag(inverse.values) %*% t(H.eigen$vectors)
  return(H.MP)
}

# matrix.Moore.Penrose2 -------------------------------------------------------
# Computes Moore-Penrose generalized inverse of symmetric matrix using 
# spectral decomposition
# By M.A.R. Ferreira
#' @export

matrix.Moore.Penrose2 <- function(H, tolp = -10) {
  H.eigen <- eigen(H)
  inverse.values <- rep(0, nrow(H))
  inverse.values[(H.eigen$values) > 10^(tolp)] <- 1 / H.eigen$values[(H.eigen$values) > 10^(tolp)]
  H.MP <- H.eigen$vectors %*% diag(inverse.values) %*% t(H.eigen$vectors)
  return(H.MP)
}

# matrix.sqrt -----------------------------------------------------------------
# Computes square root of nonnegative definite symmetric matrix using 
# spectral decomposition
#' @export

matrix.sqrt <- function(H) {
  if(nrow(H) == 1) {
    H.sqrt <- matrix(sqrt(H), nrow = 1, ncol = 1)
  } else {
    H.eigen <- eigen(H)
    H.eigen.values <- H.eigen$values    
    H.eigen.values[abs(H.eigen$values) < 10^(-20)] <- 0
    H.sqrt <- H.eigen$vectors %*% diag(sqrt(H.eigen.values)) %*% t(H.eigen$vectors)
  }  
  return(H.sqrt)
}
