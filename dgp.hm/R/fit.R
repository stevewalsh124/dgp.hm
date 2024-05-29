
# Function contents -----------------------------------------------------------
# External:
#   fit_one_layer_hm
#   fit_two_layer_hm
# Internal:
#   gibbs_one_layer_hm
#   gibbs_two_layer_hm
#   sample_w_hm
#   sample_theta_hm
#   logl_hm
#
# Updates from deepgp ---------------------------------------------------------
#   Sigma_hat is added to Sigma(X) in outer likelihood
#   No longer a nugget option (g = 1e-8 used for numerical stability)
#   No longer a prior_mean of x option
#   Matern kernel only

# fit_one_layer_hm ------------------------------------------------------------
#' @export

fit_one_layer_hm <- function(x, y, nmcmc = 10000, verb = TRUE, theta_0 = 0.1, 
                             settings = NULL, v = 2.5, Sigma_hat = NULL) {
  tic <- proc.time()[[3]]
  if (nmcmc <= 1) stop("nmcmc must be greater than 1")
  
  # Check inputs
  if (is.numeric(x)) x <- as.matrix(x)
  test <- deepgp:::check_inputs(x, y, 0) # do not need to check nugget
  settings <- deepgp:::check_settings(settings, layers = 1)
  initial <- list(theta = theta_0, tau2 = 1)
  if (!(v %in% c(0.5, 1.5, 2.5))) 
    stop("v must be one of 0.5, 1.5, or 2.5")
  if (nrow(Sigma_hat) != nrow(x))
    stop("nrow(Sigma_hat) must equal nrow(x)")
  
  out <- list(x = x, y = y, nmcmc = nmcmc, settings = settings, v = v,
              Sigma_hat = Sigma_hat)

  samples <- gibbs_one_layer_hm(x, y, nmcmc, verb, initial, settings, v, 
                                Sigma_hat)
  out <- c(out, samples)
  toc <- proc.time()[[3]]
  out$time <- unname(toc - tic)
  class(out) <- "gphm"
  return(out)
}

# fit_two_layer_hm ------------------------------------------------------------
#' @export

fit_two_layer_hm <- function (x, y, D = ifelse(is.matrix(x), ncol(x), 1), 
                              nmcmc = 10000, verb = TRUE, w_0 = NULL, 
                              theta_y_0 = 0.1, theta_w_0 = 0.1, settings = NULL, 
                              v = 2.5, Sigma_hat = NULL) {
  
  tic <- proc.time()[[3]]
  if (is.numeric(x)) x <- as.matrix(x)
  test <- deepgp:::check_inputs(x, y, 0) # do not need to check nugget
  settings <- deepgp:::check_settings(settings, layers = 2, D)
  initial <- list(w = w_0, theta_y = theta_y_0, theta_w = theta_w_0, tau2 = 1)
  initial <- deepgp:::check_initialization(initial, layers = 2, x = x, D = D)
  if (!(v %in% c(0.5, 1.5, 2.5))) 
    stop("v must be one of 0.5, 1.5, or 2.5")
  if (nrow(Sigma_hat) != nrow(x))
    stop("nrow(Sigma_hat) must equal nrow(x)")
  
  out <- list(x = x, y = y, nmcmc = nmcmc, settings = settings, v = v,
              Sigma_hat = Sigma_hat)
  
  samples <- gibbs_two_layer_hm(x, y, nmcmc, D, verb, initial, 
                                settings, v, Sigma_hat)

  out <- c(out, samples)
  toc <- proc.time()[[3]]
  out$time <- toc - tic
  class(out) <- "dgp2hm"
  return(out)
}

# gibbs_one_layer_hm ----------------------------------------------------------

gibbs_one_layer_hm <- function(x, y, nmcmc, verb, initial, settings, v, 
                               Sigma_hat) {
  dx <- deepgp::sq_dist(x)
  theta <- vector(length = nmcmc)
  theta[1] <- initial$theta
  tau2 <- vector(length = nmcmc)
  tau2[1] <- initial$tau2
  ll_store <- vector(length = nmcmc)
  ll_store[1] <- NA
  ll <- NULL
  
  for (j in 2:nmcmc) {
    
    if (verb & (j %% 500 == 0)) cat(j, "\n")
    
    # Sample lengthscale (theta), also get MLE estimate of tau2
    samp <- sample_theta_hm(y, dx, theta[j - 1], alpha = settings$alpha$theta, 
                            beta = settings$beta$theta, l = settings$l, 
                            u = settings$u, outer = TRUE, ll_prev = ll, v = v, 
                            calc_tau2 = TRUE, Sigma_hat = Sigma_hat)
    theta[j] <- samp$theta
    ll <- samp$ll
    ll_store[j] <- ll
    if (is.null(samp$tau2)) tau2[j] <- tau2[j - 1] else tau2[j] <- samp$tau2
  }
  return(list(theta = theta, tau2 = tau2, ll = ll_store))
}

# gibbs_two_layer_hm ----------------------------------------------------------

gibbs_two_layer_hm <- function(x, y, nmcmc, D, verb, initial, settings, v, 
                               Sigma_hat) {
  dx <- deepgp::sq_dist(x)
  dw <- deepgp::sq_dist(initial$w)
  theta_y <- vector(length = nmcmc)
  theta_y[1] <- initial$theta_y
  theta_w <- matrix(nrow = nmcmc, ncol = D)
  theta_w[1, ] <- initial$theta_w
  w <- list()
  w[[1]] <- initial$w
  tau2 <- vector(length = nmcmc)
  tau2[1] <- initial$tau2
  ll_store <- vector(length = nmcmc)
  ll_store[1] <- NA
  ll_outer <- NULL
  
  for (j in 2:nmcmc) {
    
    if (verb & (j %% 500 == 0)) cat(j, "\n")
    
    # Sample outer lengthscale (theta_y)
    samp <- sample_theta_hm(y, dw, theta_y[j - 1], 
                            alpha = settings$alpha$theta_y, 
                            beta = settings$beta$theta_y, l = settings$l, 
                            u = settings$u, outer = TRUE, ll_prev = ll_outer, 
                            v = v, calc_tau2 = TRUE, Sigma_hat = Sigma_hat)
    theta_y[j] <- samp$theta
    ll_outer <- samp$ll
    if (is.null(samp$tau2)) tau2[j] <- tau2[j - 1] else tau2[j] <- samp$tau2

    # Sample inner lengthscale (theta_w) - separately for each dimension
    for (i in 1:D) {
      samp <- sample_theta_hm(w[[j - 1]][, i], dx, theta_w[j - 1, i], 
                              alpha = settings$alpha$theta_w, 
                              beta = settings$beta$theta_w, l = settings$l, 
                              u = settings$u, outer = FALSE, v = v,
                              Sigma_hat = 0)
      theta_w[j, i] <- samp$theta
    }
    
    # Sample hidden Gaussian layer (w)
    samp <- sample_w_hm(y, w[[j - 1]], dw, dx, theta_y[j], theta_w[j, ], 
                        ll_prev = ll_outer, v = v, Sigma_hat = Sigma_hat)
    w[[j]] <- samp$w
    ll_outer <- samp$ll
    ll_store[j] <- ll_outer
    dw <- samp$dw
  }
  return(list(theta_y = theta_y, theta_w = theta_w, w = w, tau2 = tau2,
              ll = ll_store))
}

# sample_w_hm -----------------------------------------------------------------

sample_w_hm <- function(out_vec, w_t, w_t_dmat, in_dmat, theta_y, theta_w, 
                         ll_prev = NULL, v, Sigma_hat) {

  D <- ncol(w_t)
  if (is.null(ll_prev)) 
    ll_prev <- logl_hm(out_vec, w_t_dmat, theta_y, outer = TRUE, 
                       v = v, Sigma_hat = Sigma_hat)$logl
  for (i in 1:D) {
    
    w_prior <- mvtnorm::rmvnorm(1, sigma = deepgp:::Matern(in_dmat, 1, 
                                                           theta_w[i], 0, v))
    a <- runif(1, min = 0, max = 2 * pi)
    amin <- a - 2 * pi
    amax <- a
    ru <- runif(1, min = 0, max = 1)
    ll_threshold <- ll_prev + log(ru)
    accept <- FALSE
    count <- 0
    w_prev <- w_t[, i]
    
    while (accept == FALSE) {
      count <- count + 1
      w_t[, i] <- w_prev * cos(a) + w_prior * sin(a)
      dw <- deepgp::sq_dist(w_t)
      new_logl <- logl_hm(out_vec, dw, theta_y, outer = TRUE, 
                          v = v, Sigma_hat = Sigma_hat)$logl
      if (new_logl > ll_threshold) {
        ll_prev <- new_logl
        accept <- TRUE
      } else {
        if (a < 0) amin <- a else amax <- a
        a <- runif(1, amin, amax)
        if (count > 100) stop("reached maximum iterations of ESS")
      } # end of else statement
    } # end of while loop
  } # end of i for loop
  
  return(list(w = w_t, ll = ll_prev, dw = dw))
}

# sample_theta_hm -------------------------------------------------------------

sample_theta_hm <- function(out_vec, in_dmat, theta_t, alpha, beta, l, u, 
                            outer, ll_prev = NULL, v, calc_tau2 = FALSE, 
                            Sigma_hat) {
  
  theta_star <- runif(1, min = l * theta_t/u, max = u * theta_t/l)
  ru <- runif(1, min = 0, max = 1)
  if (is.null(ll_prev)) 
    ll_prev <- logl_hm(out_vec, in_dmat, theta_t, outer, v, 
                       Sigma_hat = Sigma_hat)$logl
  lpost_threshold <- ll_prev + dgamma(theta_t, alpha, beta, log = TRUE) + 
                        log(ru) - log(theta_t) + log(theta_star)
  ll_new <- logl_hm(out_vec, in_dmat, theta_star, outer, v, 
                    calc_tau2 = calc_tau2, Sigma_hat = Sigma_hat)
  if (ll_new$logl + dgamma(theta_star, alpha, beta, log = TRUE) > lpost_threshold) {
    return(list(theta = theta_star, ll = ll_new$logl, tau2 = ll_new$tau2))
  } else return(list(theta = theta_t, ll = ll_prev, tau2 = NULL))
}

# logl_hm ---------------------------------------------------------------------

logl_hm <- function(out_vec, in_dmat, theta, outer = FALSE, v, 
                    calc_tau2 = FALSE, Sigma_hat) {
  
  n <- length(out_vec)
  K <- deepgp:::Matern(in_dmat, 1, theta, 1e-8, v) + Sigma_hat
  id <- deepgp:::invdet(K)
  quadterm <- t(out_vec) %*% id$Mi %*% out_vec
  if (outer) {
    logl <- (-n * 0.5) * log(quadterm) - 0.5 * id$ldet
  } else logl <- (-0.5) * id$ldet - 0.5 * quadterm
  if (calc_tau2) {
    tau2 <- c(quadterm) / n
  } else tau2 <- NULL
  
  return(list(logl = c(logl), tau2 = tau2))
}
