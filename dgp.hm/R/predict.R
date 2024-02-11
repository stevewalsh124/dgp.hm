
# Function contents -----------------------------------------------------------
# External:
#   est_true.gpSW
#   est_true.dgp2SW

est_true.gpSW <- function(fit, tolpower = -10) {
  
  n <- nrow(fit$x)
  dx <- deepgp::sq_dist(fit$x)
  
  S_e <- fit$Sigma_hat
  S_ei <- moore_penrose(S_e, tolpower)
  S_ei <- (S_ei + t(S_ei)) / 2
  
  Cs <- matrix(NA, n^2, fit$nmcmc)
  Ms <- matrix(NA, n, fit$nmcmc)
  Ss <- St <- Sw <- Sx <- matrix(NA, n, fit$nmcmc)
  
  for (i in 1:fit$nmcmc) {
    
    if (i %% 500 == 0) cat(i, "\n")
    
    S <- deepgp:::Matern(dx, 1, fit$theta[i], 1e-8, fit$v)
    S_si <- solve(S)
    S_si <- (S_si + t(S_si)) / 2
    
    C <- Cs[, i] <- solve(S_si + S_ei)
    C <- (C + t(C)) / 2
    
    M <- Ms[, i] <- C %*% S_ei %*% fit$y
    
    Ss[, ((i - 1) + 1):i] <- t(mvtnorm::rmvnorm(n = 1, mean = M, sigma = C, method = "eigen"))
  }
  
  m <- rowMeans(Ss)
  ub <- apply(Ss, 1, quantile, p = 0.975)
  lb <- apply(Ss, 1, quantile, p = 0.025)
  ubb <- apply(Ss, 1, quantile, p = 0.995)
  lbb <- apply(Ss, 1, quantile, p = 0.005)
  
  return(append(fit, list(m=m, ub=ub, lb=lb, ubb=ubb, lbb=lbb, Ms=Ms, Ss=Ss, Cs=Cs)))
  
}