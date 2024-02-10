
# Function contents -----------------------------------------------------------
# External:
#   bohman
#   matrix.Moore.Penrose
#   matix.Moore.Penrose2
#   matrix.sqrt

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
