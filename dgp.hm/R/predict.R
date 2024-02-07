
# Function contents -----------------------------------------------------------
# Internal:
#   krig_SW

# krig_SW ---------------------------------------------------------------------
# Updates from deepgp:
# Remove fmin option and prior mean options
# Add Sigma_hat to Sigma(X) prior to inverting
# Add precs_pred to diagonal of predictive posterior variance

krig_SW <- function (y, dx, d_new = NULL, d_cross = NULL, theta, g, tau2 = 1, 
                     s2 = FALSE, sigma = FALSE, v = 2.5, Sigma_hat, 
                     precs_pred = rep(Inf, length(y))) {
  out <- list()
  n <- length(y)
  if (v == 999) {
    C <- deepgp:::Exp2(dx, 1, theta, g) + Sigma_hat / tau2 
    C_cross <- deepgp:::Exp2(d_cross, 1, theta, 0)
  } else {
    C <- deepgp:::Matern(dx, 1, theta, 0, v) + Sigma_hat / tau2
    C_cross <- deepgp:::Matern(d_cross, 1, theta, 0, v)
  }
  C_inv <- deepgp:::invdet(C)$Mi
  out$mean <- C_cross %*% C_inv %*% y
  
  if (s2) {
    C_new <- rep(1 + g, times = nrow(d_new)) + (1 / precs_pred) / tau2
      # ANNIE CHANGED ABOVE - added (1 / precs_pred) / tau2
    out$s2 <- tau2 * (C_new - deepgp:::diag_quad_mat(C_cross, C_inv))
  }
  
  if (sigma) {
    quadterm <- C_cross %*% C_inv %*% t(C_cross)
    if (v == 999) {
      C_new <- deepgp:::Exp2(d_new, 1, theta, g) + diag(1 / precs_pred) / tau2 # rm diag for E(Y(X))
    } else C_new <- deepgp:::Matern(d_new, 1, theta, g, v) + diag(1 / precs_pred) / tau2 # rm diag for E(Y|X)
    out$sigma <- tau2 * (C_new - quadterm)
  }
  return(out)
}
