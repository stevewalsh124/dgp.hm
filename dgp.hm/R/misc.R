
# Function contents -----------------------------------------------------------
# External:
#   opt_matern
#   scrP
# Internal:
#   nl_matern
#   moore_penrose

# Optimize Matern -------------------------------------------------------------
#' @export

opt_matern <- function(dx, y, sdd, init = c(0.1, 0.1), lower = c(-5, -5),
                       upper = c(4, 4)) { 
  out <- optim(init, nl_matern, method = "L-BFGS-B", lower = lower,
               upper = upper, dx = dx, y = y, sdd = sdd)
  return(list(theta_hat = exp(out$par[1]), 
              tau2_hat = exp(out$par[2])))
}

# Negative log likelihood Matern ----------------------------------------------

nl_matern <- function(par, dx, y, sdd) {
  # par: c(theta, tau2)
  # dx: squared distances of inputs
  # y: vector of response
  # sdd: scaled precisions
  # UPDATED: fixed kappa = 2.5 and g = 1e-8, used deepgp instead of geo
  
  theta <- exp(par[1])
  tau2 <- exp(par[2])
  n <- nrow(y)
  d <- ncol(y)
  
  K <- deepgp:::Matern(dx, tau2, theta, 1e-8, 2.5)
  id <- deepgp::invdet(K)

  ll <- 0
  for (i in 1:d) {
    yit <- y[, i] / sdd
    ll <- ll - (1/2) * (t(yit) %*% id$Mi %*% yit) - (1/2) * id$ldet
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

# moore_penrose ---------------------------------------------------------------
# Computes Moore-Penrose generalized inverse of symmetric matrix using 
# spectral decomposition
# By M.A.R. Ferreira

moore_penrose <- function(H, tolpower = -10) {
  H.eigen <- eigen(H)
  inverse.values <- rep(0, nrow(H))
  i <- (H.eigen$values > 10^(tolpower))
  inverse.values[i] <- 1 / H.eigen$values[i]
  H.MP <- H.eigen$vectors %*% diag(inverse.values) %*% t(H.eigen$vectors)
  return(H.MP)
}
