
# Function contents -----------------------------------------------------------
# External:
#   opt_matern
#   scrP
# Internal:
#   nl_matern
#   moore_penrose

# Imported Functions ----------------------------------------------------------
#' @importFrom mvtnorm rmvnorm
#' @importFrom deepgp sq_dist
#' @importFrom stats quantile dgamma optim runif
#' @importFrom grDevices heat.colors
#' @importFrom graphics lines par

# Optimize Matern -------------------------------------------------------------
#' @export

opt_matern <- function (dx, y, sdd, init = c(0.1, 0.1), n_multi = 100) {
  out <- optim(init, nl_matern, dx = dx, y = y, sdd = sdd)
  best_par <- out$par
  best_ll <- out$value
  inits <- matrix(NA,n_multi,2)
  inits[,1] <- runif(n_multi, min=init[1]-10, max=init[1]+10)
  inits[,2] <- runif(n_multi, min=init[2]-10, max=init[2]+10)
  
  for (i in 1:n_multi) {
    out <- optim(inits[i,], nl_matern, dx = dx, y = y, sdd = sdd)
    if(out$value < best_ll){
      best_par <- out$par
      best_ll <- out$value
    }
    
  }
  return(list(theta_hat = exp(best_par[2] - best_par[1]), 
              tau2_hat = exp(best_par[2])))
}

# Negative log likelihood Matern ----------------------------------------------

nl_matern <- function(par, dx, y, sdd) {
  # par: c(theta, tau2)
  # dx: squared distances of inputs
  # y: vector of response
  # sdd: scaled precisions
  # fix kappa = 2.5 and g = 1e-8
  
  # optimize (theta1, theta2) in lieu of (tau2, theta); see Walsh et. al 2023
  # theta1 = log(tau2/theta); theta2 = log(tau2)
  theta1 <- par[1]
  theta2 <- par[2]
  
  tau2 <- exp(theta2)
  theta <- exp(theta2-theta1)
  
  n <- nrow(y)
  d <- ncol(y)
  
  K <- deepgp:::Matern(dx, tau2, theta, 1e-8, 2.5)
  id <- deepgp:::invdet(K)

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

# monowarp_ref ----------------------------------------------------------------

monowarp_ref <- function(x, xg, wg, index) {
  # x: matrix of input locations for returned w
  # xg: matrix of grid locations
  # wg: UNWARPED w values at xg locations
  # Returns warped values of w at x locations
  if (!is.matrix(x)) x <- matrix(x, ncol = 1)
  if (!is.matrix(xg)) xg <- matrix(xg, ncol = 1)
  if (!is.matrix(wg)) wg <- matrix(wg, ncol = 1)
  if (!is.matrix(index)) index <- matrix(index, ncol = 1)
  
  w <- matrix(nrow = nrow(x), ncol = ncol(wg))
  for (i in 1:ncol(wg)) {
    r_orig <- range(wg[, i])
    wg[, i] <- exp(wg[, i])
    wg[, i] <- cumsum(wg[, i])
    r_new <- range(wg[, i])
    wg[, i] <- (wg[, i] - r_new[1])/(r_new[2] - r_new[1]) # first scale to [0, 1]
    wg[, i] <- wg[, i]*(r_orig[2] - r_orig[1]) + r_orig[1] # then scale to original range
    w[, i] <- deepgp:::fo_approx(xg[, i], wg[, i], x[, i], index[, i])
  }
  return(w)
}
