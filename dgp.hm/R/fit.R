
# Function contents -----------------------------------------------------------
# External:

# fit_one_layer_SW ------------------------------------------------------------
# Updates from deepgp:
# No vecchia option
# Optional Sigma_hat argument
#' @export

fit_one_layer_SW <- function(x, y, nmcmc = 10000, verb = TRUE, g_0 = 0.01, 
                             theta_0 = 0.1, true_g = NULL, settings = NULL, 
                             cov = c("matern", "exp2"), v = 2.5, 
                             Sigma_hat = NULL) {
  tic <- proc.time()[[3]]
  cov <- match.arg(cov)
  if (cov == "exp2") v <- 999 # solely used as an indicator
  if (nmcmc <= 1) stop("nmcmc must be greater than 1")
  
  # Check inputs
  if (is.numeric(x)) x <- as.matrix(x)
  test <- deepgp:::check_inputs(x, y, true_g)
  settings <- deepgp:::check_settings(settings, layers = 1)
  initial <- list(theta = theta_0, g = g_0, tau2 = 1)
  if (cov == "matern") 
    if (!(v %in% c(0.5, 1.5, 2.5))) 
      stop("v must be one of 0.5, 1.5, or 2.5")
  if (!is.null(Sigma_hat) & (length(Sigma_hat) != nrow(x)^2))
    stop("length of Sigma_hat must be the (length/nrow of x)^2")
  
  # Create output object
  out <- list(x = x, y = y, nmcmc = nmcmc, settings = settings, v = v)
  if (!is.null(Sigma_hat)) out$Sigma_hat <- Sigma_hat
  
  # Conduct MCMC
  if (is.null(Sigma_hat)) {
    samples <- deepgp:::gibbs_one_layer(x, y, nmcmc, verb, initial, true_g, 
                                        settings, v)
  } else {
    samples <- gibbs_one_layer_SW(x, y, nmcmc, verb, initial, true_g, 
                                  settings, v, Sigma_hat)
  }
  out <- c(out, samples)
  toc <- proc.time()[[3]]
  out$time <- unname(toc - tic)
  class(out) <- "gp"
  return(out)
}

# gibbs_one_layer_SW ----------------------------------------------------------
# Updates from deepgp:
# Addition of Sigma_hat to Gaussian likelihood

gibbs_one_layer_SW <- function (x, y, nmcmc, verb, initial, true_g, settings, 
                                v, Sigma_hat) {
  dx <- deepgp::sq_dist(x)
  g <- vector(length = nmcmc)
  if (is.null(true_g)) g[1] <- initial$g else g[1] <- true_g
  theta <- vector(length = nmcmc)
  theta[1] <- initial$theta
  tau2 <- vector(length = nmcmc)
  tau2[1] <- initial$tau2
  ll_store <- vector(length = nmcmc)
  ll_store[1] <- NA
  ll <- NULL
  
  for (j in 2:nmcmc) {
    if (verb & (j %% 500 == 0)) cat(j, "\n")
    
    # Sample nugget (g)
    if (is.null(true_g)) {
      samp <- sample_g_SW(y, dx, g[j - 1], theta[j - 1], 
                          alpha = settings$alpha$g, beta = settings$beta$g, 
                          l = settings$l, u = settings$u, ll_prev = ll, 
                          v = v, Sigma_hat = Sigma_hat)
      g[j] <- samp$g
      ll <- samp$ll
    } else g[j] <- true_g
    
    # Sample lengthscale (theta), also get MLE estimate of tau2
    samp <- sample_theta_SW(y, dx, g[j], theta[j - 1], alpha = settings$alpha$theta, 
                            beta = settings$beta$theta, l = settings$l, u = settings$u, 
                            outer = TRUE, ll_prev = ll, v = v, 
                            tau2 = TRUE, Sigma_hat = Sigma_hat)
    theta[j] <- samp$theta
    ll <- samp$ll
    ll_store[j] <- ll
    if (is.null(samp$tau2)) tau2[j] <- tau2[j - 1] else tau2[j] <- samp$tau2
  }
  return(list(g = g, theta = theta, tau2 = tau2, ll = ll_store))
}

# sample_g_SW -----------------------------------------------------------------
# Updates from deepgp:
# Use new log likelihood with Sigma_hat

sample_g_SW <- function (out_vec, in_dmat, g_t, theta, alpha, beta, l, u, 
                         ll_prev = NULL, v, Sigma_hat) {
  g_star <- runif(1, min = l * g_t/u, max = u * g_t/l)
  ru <- runif(1, min = 0, max = 1)
  if (is.null(ll_prev)) 
    ll_prev <- logl_SW(out_vec, in_dmat, g_t, theta, outer = TRUE, 
                       v, Sigma_hat = Sigma_hat)$logl
  lpost_threshold <- ll_prev + dgamma(g_t - deepgp:::eps, alpha, beta, log = TRUE) + 
      log(ru) - log(g_t) + log(g_star)
  ll_new <- logl_SW(out_vec, in_dmat, g_star, theta, outer = TRUE, 
                    v, Sigma_hat = Sigma_hat)$logl
  new <- ll_new + dgamma(g_star - deepgp:::eps, alpha, beta, log = TRUE)
  if (new > lpost_threshold) {
    return(list(g = g_star, ll = ll_new))
  } else return(list(g = g_t, ll = ll_prev))
}

# sample_theta_SW -------------------------------------------------------------
# Updates from deepgp: 
# Use new log likelihood with Sigma_hat

sample_theta_SW <- function (out_vec, in_dmat, g, theta_t, alpha, beta, l, u, 
                             outer, ll_prev = NULL, v, tau2 = FALSE, Sigma_hat,
                             prior_mean = 0) {
  theta_star <- runif(1, min = l * theta_t/u, max = u * theta_t/l)
  ru <- runif(1, min = 0, max = 1)
  if (is.null(ll_prev)) 
    ll_prev <- logl_SW(out_vec, in_dmat, g, theta_t, outer, v, mu = prior_mean, 
                       Sigma_hat = Sigma_hat)$logl
  lpost_threshold <- ll_prev + dgamma(theta_t, alpha, beta, log = TRUE) + 
      log(ru) - log(theta_t) + log(theta_star)
  ll_new <- logl_SW(out_vec, in_dmat, g, theta_star, outer, v, mu = prior_mean,
                    tau2 = tau2, Sigma_hat = Sigma_hat)
  if (ll_new$logl + dgamma(theta_star, alpha, beta, log = TRUE) > lpost_threshold) {
    return(list(theta = theta_star, ll = ll_new$logl, tau2 = ll_new$tau2))
  } else return(list(theta = theta_t, ll = ll_prev, tau2 = NULL))
}

# sample_w_SW -----------------------------------------------------------------
# Updates from deepgp: 
# Use new log likelihood with Sigma_hat (for acceptance, not proposal)

sample_w_SW <- function (out_vec, w_t, w_t_dmat, in_dmat, g, theta_y, theta_w, 
                         ll_prev = NULL, v, prior_mean, Sigma_hat) {
  # ANNIE CHANGED - remove cov as input
  # Can we also remove prior_mean non-zero option?
  D <- ncol(w_t)
  if (is.null(ll_prev)) 
    ll_prev <- logl_SW(out_vec, w_t_dmat, g, theta_y, outer = TRUE, 
                       v = v, Sigma_hat = Sigma_hat)$logl
  for (i in 1:D) {
    if (v == 999) {
      w_prior <- mvtnorm::rmvnorm(1, mean = prior_mean[, i], 
                                  sigma = deepgp:::Exp2(in_dmat, 1, theta_w[i], 0))
    } else w_prior <- mvtnorm::rmvnorm(1, mean = prior_mean[, i], 
                                       sigma = deepgp:::Matern(in_dmat, 1, theta_w[i], 0, v)) 
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
      new_logl <- logl_SW(out_vec, dw, g, theta_y, outer = TRUE, 
                          v = v, Sigma_hat = Sigma_hat)$logl
      if (new_logl > ll_threshold) {
        ll_prev <- new_logl
        accept <- TRUE
      } else {
        if (a < 0) {
          amin <- a
        } else {
          amax <- a
        }
        a <- runif(1, amin, amax)
        if (count > 100) stop("reached maximum iterations of ESS")
      } # end of else statement
    } # end of while loop
  } # end of i for loop
  return(list(w = w_t, ll = ll_prev, dw = dw))
}

# logl_SW ---------------------------------------------------------------------
# Updates from deepgp: 
# Add Sigma_hat to Sigma(X)

logl_SW <- function (out_vec, in_dmat, g, theta, outer = FALSE, v, 
                     tau2 = FALSE, mu = 0, Sigma_hat) {
  n <- length(out_vec)
  if (v == 999) {
    K <- deepgp:::Exp2(in_dmat, 1, theta, g) + Sigma_hat
  } else K <- deepgp:::Matern(in_dmat, 1, theta, g, v) + Sigma_hat
  id <- deepgp:::invdet(K)
  quadterm <- t(out_vec - mu) %*% id$Mi %*% (out_vec - mu)
  if (outer) {
    logl <- (-n * 0.5) * log(quadterm) - 0.5 * id$ldet
  } else logl <- (-0.5) * id$ldet - 0.5 * quadterm
  if (tau2) {
    tau2 <- c(quadterm) / n
  } else tau2 <- NULL
  
  return(list(logl = c(logl), tau2 = tau2))
}
