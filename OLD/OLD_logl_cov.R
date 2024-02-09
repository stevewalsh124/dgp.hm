# modify package functions to allow true_g to be a vector

library(deepgp)
library(Rcpp)
library(parallel)
library(doParallel)

sourceCpp("~/dgptc/dgp.hm/src/cov_SW.cpp")

##########################################
# Full likelihood function modifications #
##########################################

# link the functions from the package
check_inputs <- deepgp:::check_inputs
check_settings <- deepgp:::check_settings
check_initialization <- deepgp:::check_initialization
gibbs_two_layer <- deepgp:::gibbs_two_layer
gibbs_two_layer_vec <- deepgp:::gibbs_two_layer_vec
MaternFun <- deepgp:::Matern#Fun
Exp2Fun <- deepgp:::Exp2#Fun
invdet <- function (M){
  n <- nrow(M)
  out <- .C("inv_det_R", n = as.integer(n), M = as.double(M), 
            Mi = as.double(diag(n)), ldet = double(1))
  if(is.nan(out$ldet)) out$ldet <- c(determinant(M, logarithm = TRUE)$modulus)
  return(list(Mi = matrix(out$Mi, ncol = n), ldet = out$ldet))
}
sample_theta <- deepgp:::sample_theta
sample_theta_vec <- deepgp:::sample_theta_vec
eps <- sqrt(.Machine$double.eps)
create_approx <- deepgp:::create_approx
rand_mvn_vec <- deepgp:::rand_mvn_vec
update_obs_in_approx <- deepgp:::update_obs_in_approx


fit_two_layer_SW <- function (x, y, D = ifelse(is.matrix(x), ncol(x), 1), nmcmc = 10000, 
                              verb = TRUE, w_0 = NULL, g_0 = 0.01, theta_y_0 = 0.1, theta_w_0 = 0.1, 
                              true_g = NULL, settings = NULL, cov = c("matern", "exp2"), 
                              v = 2.5, vecchia = FALSE, m = min(25, length(y) - 1), Sigma_hat = NULL,
                              pmx = FALSE) {
  
  tic <- proc.time()[3]
  cov <- match.arg(cov)
  if (vecchia & cov == "exp2") {
    message("vecchia = TRUE requires matern covariance, proceeding with cov = 'matern'")
    cov <- "matern"
  }
  if (is.numeric(x)) 
    x <- as.matrix(x)
  test <- check_inputs(x, y, true_g)
  settings <- check_settings(settings, layers = 2, D, length(y))
  initial <- list(w = w_0, theta_y = theta_y_0, theta_w = theta_w_0, 
                  g = g_0, tau2 = 1)
  initial <- check_initialization(initial, layers = 2, x = x, 
                                  D = D, vecchia = vecchia, v = v, m = m)
  if (m >= length(y)) 
    stop("m must be less than the length of y")
  if (cov == "matern") 
    if (!(v %in% c(0.5, 1.5, 2.5))) 
      stop("v must be one of 0.5, 1.5, or 2.5")
  
  ### NEW
  if (pmx) { 
    if (ncol(x) != D) stop("for pmx, D must equal dimension of x")
    settings$w_prior_mean <- x
  } else settings$w_prior_mean <- matrix(0, nrow = nrow(x), ncol = ncol(x))
  
  out <- list(x = x, y = y, nmcmc = nmcmc, settings = settings, cov = cov)
  
  if (cov == "matern") 
    out$v <- v
  if (vecchia) 
    out$m <- m
  if (vecchia) {
    if(is.null(Sigma_hat)){
      ### NEW
      if (pmx) message("non-zero w_prior_mean is NOT implemented for this setting")
      samples <- gibbs_two_layer_vec(x, y, nmcmc, D, verb, initial, 
                                     true_g, settings, v, m)
    } else {
      if(length(Sigma_hat)==nrow(x)^2){
        samples <- gibbs_two_layer_vec_SW(x, y, nmcmc, D, verb, initial, 
                                          true_g, settings, v, m, Sigma_hat = Sigma_hat)
      } else {print(stop("length of Sigma_hat must be (length/nrow of x)^2"))}
    }
  } else {
    # check if length(Sigma_hat) is NULL, or nrow(as.matrix(x))
    if(is.null(Sigma_hat)){
      ### NEW
      if (pmx) message("non-zero w_prior_mean is NOT implemented for this setting")
      samples <- gibbs_two_layer(x, y, nmcmc, D, verb, initial, 
                                 true_g, settings, cov, v)
    } else {
      if(length(Sigma_hat)==nrow(x)^2){
        samples <- gibbs_two_layer_SW(x, y, nmcmc, D, verb, initial, 
                                      true_g, settings, cov, v, Sigma_hat)
      } else {print(stop("length of Sigma_hat must be the (length/nrow of x)^2"))}
    }
  }
  out <- c(out, samples)
  toc <- proc.time()[3]
  out$time <- toc - tic
  if(!is.null(Sigma_hat)) 
    out$Sigma_hat <- Sigma_hat
  if (vecchia) 
    class(out) <- "dgp2vec"
  else class(out) <- "dgp2"
  return(out)
}


# use this version when g is a vector
gibbs_two_layer_SW <- function (x, y, nmcmc, D, verb, initial, true_g, settings, cov, v, Sigma_hat){
  print("using SW")
  dx <- sq_dist(x)
  dw <- sq_dist(initial$w)
  g <- vector(length = nmcmc)
  if (is.null(true_g)) 
    g[1] <- initial$g
  else g[1] <- true_g
  theta_y <- vector(length = nmcmc)
  theta_y[1] <- initial$theta_y
  theta_w <- matrix(nrow = nmcmc, ncol = D)
  theta_w[1, ] <- initial$theta_w
  w <- list()
  w[[1]] <- initial$w
  tau2 <- vector(length = nmcmc)
  tau2[1] <- initial$tau2
  ll_outer <- NULL
  for (j in 2:nmcmc) {
    if (verb) 
      if (j%%1000 == 0) 
        cat(j, "\n")
    if (is.null(true_g)) {
      samp <- sample_g_SW(y, dw, g[j - 1], theta_y[j - 1], 
                          alpha = settings$alpha$g, beta = settings$beta$g, 
                          l = settings$l, u = settings$u, ll_prev = ll_outer, 
                          v = v, cov = cov, Sigma_hat = Sigma_hat)
      g[j] <- samp$g
      ll_outer <- samp$ll
      # NOTE - outer layer unaffected by pmx
    } else {g[j] <- true_g}
    samp <- sample_theta_SW(y, dw, g[j], theta_y[j - 1], alpha = settings$alpha$theta_y, 
                            beta = settings$beta$theta_y, l = settings$l, u = settings$u, 
                            outer = TRUE, ll_prev = ll_outer, v = v, cov = cov, 
                            tau2 = TRUE, Sigma_hat = Sigma_hat)
    theta_y[j] <- samp$theta
    ll_outer <- samp$ll
    if (is.null(samp$tau2)) 
      tau2[j] <- tau2[j - 1]
    else tau2[j] <- samp$tau2
    # NOTE - outer layer unaffected by pmx
    
    for (i in 1:D) {
      ### NEW - sample_theta_new function accepts prior mean as argument
      ### NEW - calls logl_new which accepts prior mean as argument
      samp <- sample_theta_new(w[[j - 1]][, i], dx, g = eps, 
                           theta_w[j - 1, i], alpha = settings$alpha$theta_w, 
                           beta = settings$beta$theta_w, l = settings$l, 
                           u = settings$u, outer = FALSE, v = v,
                           prior_mean = settings$w_prior_mean[, i])
      theta_w[j, i] <- samp$theta
    }
    
    samp <- sample_w_SW(y, w[[j - 1]], dw, dx, g[j], theta_y[j], 
                        theta_w[j, ], ll_prev = ll_outer, v = v, cov = cov, 
                        prior_mean = settings$w_prior_mean, Sigma_hat = Sigma_hat)
    w[[j]] <- samp$w
    ll_outer <- samp$ll
    dw <- samp$dw
  }
  return(list(g = g, theta_y = theta_y, theta_w = theta_w, 
              w = w, tau2 = tau2))
}


