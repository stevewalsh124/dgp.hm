
# Function contents -----------------------------------------------------------
# External:
#   est_true (applies to both gpSW and dgp2SW classes)

# est_true --------------------------------------------------------------------
#' @export 

est_true <- function(fit, tolpower = -10) {
  # this function finds the posterior mean of the true spectrum
  # given the observed weighted average and precision information
  n <- nrow(fit$x)
  if (!(class(fit) %in% c("gpSW", "dgp2SW")))
    stop("only applicable to gpSW or dgp2SW classes")
  
  deep <- (class(fit) == "dgp2SW")
  
  # get squared distance matrix
  if (!deep) d2 <- deepgp::sq_dist(fit$x)
  
  # From Section 3.3.2: Bayesian UQ Algorithm (Walsh dissertation 2023)
  # Get Sigma_epsilon (covariance matrix for correlated errors) and its inverse
  # Cov mtx of [Ybar]
  S_e <- fit$Sigma_hat
  S_ei <- moore_penrose(S_e, tolpower)
  S_ei <- (S_ei + t(S_ei)) / 2
  
  # Store the Covariance matrices/Means/Simulated curves for each mcmc iter
  Cs <- matrix(NA, n^2, fit$nmcmc)
  Ms <- matrix(NA, n, fit$nmcmc)
  Ss <- St <- Sw <- Sx <- matrix(NA, n, fit$nmcmc)
  
  # For each mcmc iteration...
  for (i in 1:fit$nmcmc) {
    
    if (i %% 500 == 0) cat(i, "\n")
    
    # construct the covariance matrix for the true spectrum (Sigma_S)
    # cov mtx of [S]; based on a matern covariance function
    if (deep) {
      d2 <- deepgp:::sq_dist(fit$w[[i]])
      theta <- fit$theta_y[i]
    } else theta <- fit$theta[i]
    S <- deepgp:::Matern(d2, 1, theta, 1e-8, fit$v)
    
    # Find inverse of Sigma_S, and ensure it's symmetric
    S_si <- solve(S)
    S_si <- (S_si + t(S_si)) / 2
    
    # Find the cov matrix for the true spectrum, given the observed avg 
    # cov mtx of [S | Ybar]
    C <- Cs[, i] <- solve(S_si + S_ei)
    C <- (C + t(C)) / 2
    
    # The mean of [S|ybar] depends on C
    M <- Ms[, i] <- C %*% S_ei %*% fit$y
    
    # For these values of S and C (which vary by mcmc iter),
    # Simulate a realization from the distbn [S|Ybar] ~ N(M,C)
    Ss[, ((i - 1) + 1):i] <- t(mvtnorm::rmvnorm(n = 1, mean = M, sigma = C, 
                                                method = "eigen"))
  }
  
  # Get the posterior mean of [S|Ybar], as well as credible intervals
  fit$m <- rowMeans(Ss)
  fit$ub <- apply(Ss, 1, quantile, p = 0.975)
  fit$lb <- apply(Ss, 1, quantile, p = 0.025)
  fit$ubb <- apply(Ss, 1, quantile, p = 0.995)
  fit$lbb <- apply(Ss, 1, quantile, p = 0.005)
  # fit$Ms <- Ms; fit$Ss <- Ss; fit$Cs <- Cs
  
  return(fit)
  
}
