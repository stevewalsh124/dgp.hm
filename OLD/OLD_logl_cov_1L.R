# modify package functions to allow true_g to be a vector

#################################
# Prediction code modifications #
#################################

clean_prediction <- deepgp:::clean_prediction

predict.dgp2_SW <- function (object, x_new, lite = TRUE, store_latent = FALSE, mean_map = TRUE, 
                             EI = FALSE, cores = detectCores() - 1, precs_pred = NULL, ...){
  tic <- proc.time()[3]
  # print(object$cov); print(object$v); print(object$nmcmc)
  # object <- clean_prediction(object)
  # print(object$cov); print(object$v); print(object$nmcmc)
  if(object$cov == "exp2"){object$v <- 999}
  if (is.numeric(x_new)) 
    x_new <- as.matrix(x_new)
  object$x_new <- x_new
  n_new <- nrow(object$x_new)
  D <- ncol(object$w[[1]])
  dx <- sq_dist(object$x)
  d_new <- sq_dist(object$x_new)
  d_cross <- sq_dist(object$x_new, object$x)
  iters <- 1:object$nmcmc
  if (cores == 1) {
    chunks <- list(iters)
  }
  else chunks <- split(iters, sort(cut(iters, cores, labels = FALSE)))
  if (cores > detectCores()) 
    warning("cores is greater than available nodes")
  cl <- makeCluster(cores)
  clusterExport(cl, c("krig_SW", "MaternFun", "eps", "invdet", "sq_dist", "Exp2Fun"))
  registerDoParallel(cl)
  thread <- NULL
  result <- foreach(thread = 1:cores) %dopar% {
    out <- list()
    if (store_latent) 
      out$w_new <- list()
    out$mu_t <- matrix(nrow = n_new, ncol = length(chunks[[thread]]))
    if (lite) {
      out$s2_sum <- rep(0, times = n_new)
    }
    else out$sigma_sum <- matrix(0, nrow = n_new, ncol = n_new)
    if (EI) 
      out$ei_sum <- rep(0, times = n_new)
    j <- 1
    for (t in chunks[[thread]]) {
      w_t <- object$w[[t]]
      w_new <- matrix(nrow = n_new, ncol = D)
      for (i in 1:D) {
        if (mean_map) {
          k <- krig_SW(w_t[, i], dx, NULL, d_cross, object$theta_w[t, i], 
                       g = eps, v = object$v, Sigma_hat = object$Sigma_hat, precs_pred = precs_pred)
          w_new[, i] <- k$mean
        }
        else {
          k <- krig_SW(w_t[, i], dx, d_new, d_cross, object$theta_w[t, i], 
                       g = eps, sigma = TRUE, v = object$v, Sigma_hat = object$Sigma_hat, precs_pred = precs_pred)
          w_new[, i] <- mvtnorm::rmvnorm(1, k$mean, k$sigma)
        }
      }
      if (store_latent) 
        out$w_new[[j]] <- w_new
      k <- krig_SW(object$y, sq_dist(w_t), sq_dist(w_new), 
                   sq_dist(w_new, w_t), object$theta_y[t], object$g[t], 
                   object$tau2[t], s2 = lite, sigma = !lite, f_min = EI, 
                   v = object$v, Sigma_hat = object$Sigma_hat, precs_pred = precs_pred)
      out$mu_t[, j] <- k$mean
      if (lite) {
        out$s2_sum <- out$s2_sum + k$s2
      }
      else out$sigma_sum <- out$sigma_sum + k$sigma
      if (EI) {
        if (lite) {
          sig2 <- k$s2 - (object$tau2[t] * object$g[t])
        }
        else sig2 <- diag(k$sigma) - (object$tau2[t] * object$g[t])
        out$ei_sum <- out$ei_sum + exp_improv(k$mean, sig2, k$f_min)
      }
      j <- j + 1
    }
    return(out)
  }
  stopCluster(cl)
  mu_t <- do.call(cbind, lapply(result, with, eval(parse(text = "mu_t"))))
  if (lite) {
    s2_sum <- Reduce("+", lapply(result, with, eval(parse(text = "s2_sum"))))
  }
  else {
    sigma_sum <- Reduce("+", lapply(result, with, eval(parse(text = "sigma_sum"))))
  }
  if (store_latent) 
    w_new <- unlist(lapply(result, with, eval(parse(text = "w_new"))), 
                    recursive = FALSE)
  if (EI) 
    ei_sum <- Reduce("+", lapply(result, with, eval(parse(text = "ei_sum"))))
  mu_cov <- cov(t(mu_t))
  object$mean <- rowMeans(mu_t)
  if (store_latent) 
    object$w_new <- w_new
  if (lite) {
    object$s2 <- s2_sum/object$nmcmc + diag(mu_cov)
    object$s2_smooth <- object$s2 - mean(object$g * object$tau2)/precs_pred
  }
  else {
    object$Sigma <- sigma_sum/object$nmcmc + mu_cov
    object$Sigma_smooth <- object$Sigma - diag(mean(object$g * object$tau2)/precs_pred, n_new)
  }
  if (EI) 
    object$EI <- ei_sum/object$nmcmc
  toc <- proc.time()[3]
  object$time <- object$time + (toc - tic)
  return(object)
}

