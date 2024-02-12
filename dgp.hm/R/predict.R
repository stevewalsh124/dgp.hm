
# Function contents -----------------------------------------------------------
# External:
#   est_true (applies to both gpSW and dgp2SW classes)

# est_true --------------------------------------------------------------------
#' @export 

est_true <- function(fit, tolpower = -10) {
  
  n <- nrow(fit$x)
  if (!(class(fit) %in% c("gpSW", "dgp2SW")))
    stop("only applicable to gpSW or dgp2SW classes")
  
  deep <- (class(fit) == "dgp2SW")
  if (!deep) d2 <- deepgp::sq_dist(fit$x)
  
  S_e <- fit$Sigma_hat
  S_ei <- moore_penrose(S_e, tolpower)
  S_ei <- (S_ei + t(S_ei)) / 2
  
  Cs <- matrix(NA, n^2, fit$nmcmc)
  Ms <- matrix(NA, n, fit$nmcmc)
  Ss <- St <- Sw <- Sx <- matrix(NA, n, fit$nmcmc)
  
  for (i in 1:fit$nmcmc) {
    
    if (i %% 500 == 0) cat(i, "\n")
    
    if (deep) {
      d2 <- deepgp:::sq_dist(fit$w[[i]])
      theta <- fit$theta_y[i]
    } else theta <- fit$theta[i]
    
    S <- deepgp:::Matern(d2, 1, theta, 1e-8, fit$v)
    S_si <- solve(S)
    S_si <- (S_si + t(S_si)) / 2
    
    C <- Cs[, i] <- solve(S_si + S_ei)
    C <- (C + t(C)) / 2
    
    M <- Ms[, i] <- C %*% S_ei %*% fit$y
    
    Ss[, ((i - 1) + 1):i] <- t(mvtnorm::rmvnorm(n = 1, mean = M, sigma = C, 
                                                method = "eigen"))
  }
  
  fit$m <- rowMeans(Ss)
  fit$ub <- apply(Ss, 1, quantile, p = 0.975)
  fit$lb <- apply(Ss, 1, quantile, p = 0.025)
  fit$ubb <- apply(Ss, 1, quantile, p = 0.995)
  fit$lbb <- apply(Ss, 1, quantile, p = 0.005)
  # fit$Ms <- Ms; fit$Ss <- Ss; fit$Cs <- Cs
  
  return(fit)
  
}
