
# Function contents -----------------------------------------------------------
# External:


# logl_SW ---------------------------------------------------------------------
# Updates from deepgp: 
# Add Sigma_hat to Sigma(X)

logl_SW <- function (out_vec, in_dmat, g, theta, outer = FALSE, v, 
                     tau2 = FALSE, Sigma_hat) {
  # ANNIE CHANGED - remove cov as input
  n <- length(out_vec)
  if (v = 999) {
    K <- deepgp:::Matern(in_dmat, 1, theta, g, v) + Sigma_hat
  } else K <- deepgp:::Exp2(in_dmat, 1, theta, g) + Sigma_hat
  id <- invdet(K)
  quadterm <- t(out_vec) %*% id$Mi %*% (out_vec)
  if (outer) {
    logl <- (-n * 0.5) * log(quadterm) - 0.5 * id$ldet
  } else logl <- (-0.5) * id$ldet - 0.5 * quadterm
  if (tau2) {
    tau2 <- 1#c(quadterm)/n
    # WHY NOT ESTIMATE TAU2?
  } else tau2 <- NULL
  
  return(list(logl = c(logl), tau2 = tau2))
}
