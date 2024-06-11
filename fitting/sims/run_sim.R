###############################################################################
# This script conducts a "bake-off" of various surrogate models under 2
# different test functions (f1 and f2) and 4 different noise settings for 
# a specified random seed.
#
# Important Arguments:
#   seed: random seed (must be set before sourcing functions.R)
#   func: 1 or 2
#   setting: 1, 2, 3, or 4 (indicates noise settings A-D from DPC paper)
#   r: number of function realizations
#
# Output:
#   data frame of MSE values for each model is either written to ".csv" in the
#   "results" folder or appended to existing file
#   
# Included Models:
#   dgp.hm: the Bayesian hierarchical DGP model 
#   deepgp: original DGP from "deepgp" package
#   hetgp: heteroskedastic GP from "hetGP" package 

# Note: only dgp.hm is actually designed for multiple response observations
###############################################################################

library(dgp.hm) # must install locally
library(deepgp)
library(hetGP)
library(mvtnorm)

seed <- 1
func <- 1
setting <- 3
r <- 5
args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))
if(!(func %in% 1:2)) stop("func should be 1 or 2")
if(!(setting %in% 1:4)) stop("setting should be 1, 2, 3, or 4")

set.seed(seed)
source("functions.R") # generates random model parameters upon sourcing

vis <- FALSE # should plots be generated

# Generate data ---------------------------------------------------------------

x <-seq(0, 4, by=0.1)
n <- length(x)
Sigma_true <- get_Sigma_true(x, n, func, setting)
if(vis) image(Sigma_true) # make sure sd plot looks right

Y <- matrix(nrow = r, ncol = n)
if(func == 1) {
  for (i in 1:r) Y[i,] <- f1(x, m1=m1, u1=u1, Sigma = Sigma_true)
} else {
  for (i in 1:r) Y[i,] <- f2(x, m2=m2, u2=u2, Sigma = Sigma_true)
}
if(vis) matplot(x, t(Y), ylab="f(x)", type="l")

y_avg <- colMeans(Y)
if(func == 1) y_true <- f1(x, m1=m1, u1=u1)
if(func == 2) y_true <- f2(x, m2=m2, u2=u2)

# x_all <- NULL
# y_all <- NULL
# for(i in 1:r) {
#   x_all <- c(x_all, x)
#  y_all <- c(y_all, Y[i, ])
# }

# dgp.hm model ----------------------------------------------------------------

# Previous option - estimate Sigma_hat (not used currently)
# if(setting %in% 1:3) {
#   var_y <- mean(apply(Y, 2, var))
#   Sigma_hat = diag(var_y, n)
# } else {
#   params_hat <- opt_matern(D, t(Y), sdd =rep(1,n))
#   Sigma_hat <- deepgp:::Matern(D, params_hat$tau2_hat, params_hat$theta_hat,
#                                g=1e-8, v=2.5)
# }

fit <- dgp.hm::fit_two_layer_hm(x, y_avg, Sigma_hat = Sigma_true, nmcmc = 20000)
if(vis) plot(fit)
fit <- dgp.hm::trim(fit, 15000, 5)
fit <- dgp.hm::est_true(fit)

if (vis) {
  matplot(x, t(Y), type="l", col="gray")
  lines(x, y_true)
  lines(x, y_avg, col="red", lty=2)
  lines(x, fit$m, col="blue")
  lines(x, fit$ub, col="blue")
  lines(x, fit$lb, col="blue")
  legend(x = "topright", legend = c("data","truth", "wt avg", "95% UQ"),
         col = c("gray","black","red","blue"), lty = c(1,1,2,1))
}
dgp_mse <- mean((fit$m - y_true)^2)

# deepgp model ----------------------------------------------------------------
# Not quite sure what to do here - use the average?  Or all the data?

fit2 <- deepgp::fit_two_layer(x, y_avg, nmcmc = 10000, 
                              true_g = mean(diag(Sigma_true)))
if(vis) plot(fit2)
fit2 <- deepgp::trim(fit2, 5000, 5)
fit2 <- predict(fit2, x)
if(vis) plot(fit2)
deepgp_mse <- mean((fit2$mean - y_true)^2)

# hetGP model -----------------------------------------------------------------
# Not sure how to pass the true variance to this model....

fit3 <- mleHetGP(matrix(x_all, ncol = 1), y_all, covtype = "Matern5_2")
pred <- predict(matrix(x, ncol = 1), object = fit3)
if (vis) {
  matplot(x, t(Y), type="l", col="gray")
  lines(x, y_true)
  lines(x, y_avg, col="red", lty=2)
  lines(x, pred$mean, col="blue")
  lines(x, pred$mean - 2*sqrt(pred$sd2+pred$nugs), col="blue")
  lines(x, pred$mean + 2*sqrt(pred$sd2+pred$nugs), col="blue")
  legend(x = "topright", legend = c("data","truth", "wt avg", "95% UQ"),
         col = c("gray","black","red","blue"), lty = c(1,1,2,1))
}
hetgp_mse <- mean((pred$mean - y_true)^2)

# Store results ---------------------------------------------------------------

filename <- paste0("results/sims_", func, "_", setting, "_", r, ".csv")
if (file.exists(filename)) {
  results <- read.csv(filename)
  results <- rbind(results, c(seed, dgp_mse, deepgp_mse, hetgp_mse))
} else {
  results <- data.frame(seed = seed, dgp = dgp_mse, deepgp = deepgp_mse,
                        hetgp = hetgp_mse)
}
write.csv(results, filename, row.names = FALSE)

