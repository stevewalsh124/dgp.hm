
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
#   one dimension only
#   Matern kernel only
#   scale tau2 is fixed at 1
#   No longer a nugget option (g = 1e-8 used for numerical stability)
#   No longer a prior_mean of x option

# fit_one_layer_hm ------------------------------------------------------------
#' @export

fit_one_layer_hm <- function(x, y, Sigma_hat, nmcmc = 10000, verb = TRUE, theta_0 = 0.1, 
                             settings = NULL, v = 2.5) {
  tic <- proc.time()[[3]]
  if (nmcmc <= 1) stop("nmcmc must be greater than 1")
  
  # Check inputs
  if (!is.matrix(x)) {
    x <- matrix(x, ncol = 1)
  } else if (ncol(x) != 1) stop("x must be one-dimensional")
  if (nrow(x) != length(y)) stop("dimensions of x and y do not match")
  settings <- check_prior_settings(settings, layers = 1)
  initial <- list(theta = theta_0)
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

fit_two_layer_hm <- function (x, y, Sigma_hat, nmcmc = 10000, verb = TRUE, 
                              w_0 = NULL, theta_y_0 = 0.1, theta_w_0 = 0.1, 
                              settings = NULL, v = 2.5, x_grid = NULL) {
  
  tic <- proc.time()[[3]]
  if (!is.matrix(x)) {
    x <- matrix(x, ncol = 1)
  } else if (ncol(x) != 1) stop("x must be one-dimensional")
  if (nrow(x) != length(y)) stop("dimensions of x and y do not match")
  settings <- check_prior_settings(settings, layers = 2)
  if (is.null(w_0)) w_0 <- x
  initial <- list(w = w_0, theta_y = theta_y_0, theta_w = theta_w_0)
  if (is.null(x_grid)) 
    x_grid <- matrix(seq(0, 1, length = 50), ncol = 1)
  if (!(v %in% c(0.5, 1.5, 2.5))) 
    stop("v must be one of 0.5, 1.5, or 2.5")
  if (nrow(Sigma_hat) != nrow(x))
    stop("nrow(Sigma_hat) must equal nrow(x)")
  
  out <- list(x = x, y = y, x_grid = x_grid, nmcmc = nmcmc, settings = settings, 
              v = v, Sigma_hat = Sigma_hat)
  
  samples <- gibbs_two_layer_hm(x, y, x_grid, nmcmc, verb, initial, 
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
  ll_store <- vector(length = nmcmc)
  ll_store[1] <- NA
  ll <- NULL
  
  for (j in 2:nmcmc) {
    if (verb & (j %% 500 == 0)) cat(j, "\n")
    samp <- sample_theta_hm(y, dx, theta[j - 1], alpha = settings$theta$alpha, 
                            beta = settings$theta$beta, l = settings$l, 
                            u = settings$u, ll_prev = ll, v = v, 
                            Sigma_hat = Sigma_hat)
    theta[j] <- samp$theta
    ll <- samp$ll
    ll_store[j] <- ll
  }
  
  return(list(theta = theta, ll = ll_store))
}

# gibbs_two_layer_hm ----------------------------------------------------------

gibbs_two_layer_hm <- function(x, y, x_grid, nmcmc, verb, initial, settings, v, 
                               Sigma_hat) {
  
  ng <- nrow(x_grid)
  dx_grid <- deepgp::sq_dist(x_grid)
  grid_index <- deepgp:::fo_approx_init(x_grid, x)
  
  theta_y <- vector(length = nmcmc)
  theta_y[1] <- initial$theta_y
  theta_w <- vector(length = nmcmc)
  theta_w[1] <- initial$theta_w
  w_grid <- matrix(nrow = ng, ncol = nmcmc)
  w_grid[, 1] <- deepgp:::fo_approx(x, initial$w, x_grid) # snap initial$w to grid
  w <- matrix(nrow = nrow(x), ncol = nmcmc)
  w[, 1] <- monowarp_ref(x, x_grid, w_grid[, 1], grid_index)
  ll_store <- vector(length = nmcmc)
  ll_store[1] <- NA
  ll_outer <- NULL
  
  dw <- deepgp::sq_dist(w[, 1])
  
  for (j in 2:nmcmc) {
    if (verb & (j %% 500 == 0)) cat(j, "\n")
    
    # Sample outer lengthscale (theta_y)
    samp <- sample_theta_hm(y, dw, theta_y[j - 1], 
                            alpha = settings$theta_y$alpha, 
                            beta = settings$theta_y$beta, l = settings$l, 
                            u = settings$u, ll_prev = ll_outer, 
                            v = v, Sigma_hat = Sigma_hat)
    theta_y[j] <- samp$theta
    ll_outer <- samp$ll

    # Sample inner lengthscale (theta_w)
    samp <- sample_theta_hm(w_grid[, j - 1], dx_grid, theta_w[j - 1], 
                            alpha = settings$theta_w$alpha, 
                            beta = settings$theta_w$beta, l = settings$l, 
                            u = settings$u, v = v, Sigma_hat = 0)
    theta_w[j] <- samp$theta
    
    # Sample hidden Gaussian layer (w)
    samp <- sample_w_hm(y, w_grid[, j - 1], dw, x, x_grid, dx_grid, grid_index,
                        theta_y[j], theta_w[j], ll_prev = ll_outer, v = v, 
                        Sigma_hat = Sigma_hat)
    w_grid[, j] <- samp$w_grid
    w[, j] <- samp$w
    dw <- samp$dw
    ll_outer <- samp$ll
    ll_store[j] <- ll_outer
  }
  
  return(list(theta_y = theta_y, theta_w = theta_w, w = w, w_grid = w_grid, 
              ll = ll_store))
}

# sample_w_hm -----------------------------------------------------------------

sample_w_hm <- function(y, w_grid, dw, x, x_grid, dx_grid, grid_index, 
                        theta_y, theta_w, ll_prev = NULL, v, Sigma_hat) {

  if (is.null(ll_prev)) 
    ll_prev <- logl_hm(y, dw, theta_y, v, Sigma_hat)
  
  w_grid_prior <- mvtnorm::rmvnorm(1, sigma = deepgp:::Matern(dx_grid, 1, theta_w, 0, v))
  a <- runif(1, min = 0, max = 2*pi)
  amin <- a - 2*pi
  amax <- a
  ru <- runif(1, min = 0, max = 1)
  ll_threshold <- ll_prev + log(ru)
  accept <- FALSE
  count <- 0

  while (accept == FALSE) {
    count <- count + 1
    w_grid_new <- w_grid*cos(a) + w_grid_prior*sin(a)
    w <- monowarp_ref(x, x_grid, w_grid_new, grid_index)
    dw <- deepgp::sq_dist(w)
    ll_new <- logl_hm(y, dw, theta_y, v, Sigma_hat)
    if (ll_new > ll_threshold) {
      accept <- TRUE
    } else {
      if (a < 0) amin <- a else amax <- a
      a <- runif(1, amin, amax)
      if (count > 100) stop("reached maximum iterations of ESS")
    } # end of else statement
  } # end of while loop

  return(list(w_grid = w_grid_new, w = w, dw = dw, ll = ll_new))
}

# sample_theta_hm -------------------------------------------------------------

sample_theta_hm <- function(y, dx, theta_t, alpha, beta, l, u, 
                            ll_prev = NULL, v, Sigma_hat) {
  
  theta_star <- runif(1, min = l * theta_t/u, max = u * theta_t/l)
  ru <- runif(1, min = 0, max = 1)
  if (is.null(ll_prev)) 
    ll_prev <- logl_hm(y, dx, theta_t, v, Sigma_hat)
  lpost_threshold <- ll_prev + dgamma(theta_t, alpha, beta, log = TRUE) + 
                        log(ru) - log(theta_t) + log(theta_star)
  ll_new <- logl_hm(y, dx, theta_star, v, Sigma_hat)
  if (ll_new + dgamma(theta_star, alpha, beta, log = TRUE) > lpost_threshold) {
    return(list(theta = theta_star, ll = ll_new))
  } else return(list(theta = theta_t, ll = ll_prev))
}

# logl_hm ---------------------------------------------------------------------

logl_hm <- function(y, dx, theta, v, Sigma_hat) {
  
  n <- length(y)
  K <- deepgp:::Matern(dx, 1, theta, 1e-8, v) + Sigma_hat
  id <- deepgp:::invdet(K)
  quadterm <- t(y) %*% id$Mi %*% y
  ll <- (-0.5)*id$ldet - 0.5*quadterm
  return(ll)
}
