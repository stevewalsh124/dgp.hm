###############################################################################
# This script conducts a "bake-off" of various surrogate models under 2
# different test functions (f1 and f2) and 4 different noise settings for 
# a specified random seed.  All models are tasked with estimating the noise.
#
# Important Arguments:
#   seed: random seed (must be set before sourcing functions.R)
#   func: 1 or 2
#   setting: 1, 2, 3, or 4 (indicates noise settings A-D from DPC paper)
#
# Output:
#   data frame of MSE and CRPS values for each model is either written to 
#   ".csv" in the "results" folder or appended to existing file
#   
# Included Models:
#   dgp: the Bayesian hierarchical DGP model 
#   dgp_r: the dgp model, but use Sigma_hat/sqrt(r)
#   deepgp: original DGP from "deepgp" package
#   hetgp: heteroskedastic GP from "hetGP" package 

# Note: only dgp.hm is actually designed for multiple response observations
###############################################################################

# TODO: bump up the reps from 5 to 20, focus on setting 4,
# ANNIE TODO: code up coverage and logs for each method, rerun for r = 5 and r = 20
# might need to make the Sigma_true a bit more complicated
# Note to Annie: adjust default priors like I did for the deepgp package

library(dgp.hm) # must install locally
library(deepgp)
library(hetGP)
library(mvtnorm)
library(scoringRules)

seed <- 1
func <- 1
setting <- 1
args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))
if(!(func %in% 1:2)) stop("func should be 1 or 2")
if(!(setting %in% 1:4)) stop("setting should be 1, 2, 3, or 4")

set.seed(seed)
source("functions.R") # generates random model parameters upon sourcing

r <- 5 # number of replicates observed
vis <- FALSE # should plots be generated

# Generate data ---------------------------------------------------------------

x <- seq(0, 4, by = 0.1)
n <- length(x)
Sigma_true <- get_Sigma_true(x, n, func, setting)
if(vis) image(Sigma_true) # make sure sd plot looks right

Y <- matrix(nrow = r, ncol = n)
if(func == 1) {
  for (i in 1:r) Y[i,] <- f1(x, m1 = m1, u1 = u1, Sigma = Sigma_true)
} else {
  for (i in 1:r) Y[i,] <- f2(x, m2 = m2, u2 = u2, Sigma = Sigma_true)
}
if(vis) matplot(x, t(Y), ylab = "f(x)", type = "l")

if(func == 1) y_true <- f1(x, m1 = m1, u1 = u1, Sigma = diag(1e-300,n,n))
if(func == 2) y_true <- f2(x, m2 = m2, u2 = u2, Sigma = diag(1e-300,n,n))

# get zero-mean version of Y
y_avg <- colMeans(Y)
Yzm <- matrix(NA, r, n)
for (i in 1:r) Yzm[i,] <- Y[i,] - y_avg

x_all <- NULL
y_all <- NULL
for(i in 1:r) {
  x_all <- c(x_all, x)
  y_all <- c(y_all, Y[i, ])
}
if(vis) plot(x_all, y_all)

# dgp.hm model ----------------------------------------------------------------

tic <- proc.time()[3]
if(setting %in% 1:3) {
  var_y <- mean(apply(Y, 2, var))
  Sigma_hat <- diag(var_y, n)
} else {
  dx <- sq_dist(x)
  params_hat <- opt_matern(dx, t(Yzm), sdd = rep(1,n), 
                           init = c(((max(x)-min(x))/10), mean(apply(Yzm, 2, var))))
  Sigma_hat <- deepgp:::Matern(dx, params_hat$tau2_hat, params_hat$theta_hat,
                               g = 1e-8, v = 2.5)
}

fit <- dgp.hm::fit_two_layer_hm(x, y_avg, Sigma_hat = Sigma_hat, 
                                nmcmc = 10000)
fit <- dgp.hm::trim(fit, 5000, 5)
if(vis) plot(fit)
fit <- dgp.hm::est_true(fit, return_all = TRUE)

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
# 1. check coverage??? 95% of the samples, do they contain the truth
# 2. log score, assuming normal distribution, estimate sd of the samples
dgp_logs <- sum(logs.numeric(drop(y_true), family = "normal", mean = fit$m, 
                             sd = apply(fit$Ss, 1, sd)))
toc <- proc.time()[3]
dgp_time <- toc - tic

# dgp_r model --------------------------------------------------------------

tic <- proc.time()[3]
fit1a <- dgp.hm::fit_two_layer_hm(x, y_avg, Sigma_hat = Sigma_hat/sqrt(r), nmcmc = 10000)
if(vis) plot(fit1a)
fit1a <- dgp.hm::trim(fit1a, 5000, 5)
fit1a <- dgp.hm::est_true(fit1a)

if (vis) {
  matplot(x, t(Y), type="l", col="gray")
  lines(x, y_true)
  lines(x, y_avg, col="red", lty=2)
  lines(x, fit1a$m, col="blue")
  lines(x, fit1a$ub, col="blue")
  lines(x, fit1a$lb, col="blue")
  legend(x = "topright", legend = c("data","truth", "wt avg", "95% UQ"),
         col = c("gray","black","red","blue"), lty = c(1,1,2,1))
}
dgp_r_mse <- mean((fit1a$m - y_true)^2)
toc <- proc.time()[3]
time_r <- toc - tic

# deepgp model ----------------------------------------------------------------

tic <- proc.time()[3]
fit2 <- deepgp::fit_two_layer(x_all, y_all, nmcmc = 10000, vecchia = T)
if(vis) plot(fit2)
fit2 <- deepgp::trim(fit2, 5000, 5)
fit2 <- predict(fit2, x)
if(vis) plot(fit2)
deepgp_mse <- mean((fit2$mean - y_true)^2)
toc <- proc.time()[3]
time2 <- toc - tic

# hetGP model -----------------------------------------------------------------

tic <- proc.time()[3]
fit3 <- mleHetGP(matrix(x_all, ncol = 1), y_all, covtype = "Matern5_2")
pred <- predict(matrix(x, ncol = 1), object = fit3)
if (vis) {
  matplot(x, t(Y), type="l", col="gray")
  lines(x, y_true)
  lines(x, y_avg, col="red", lty=2)
  lines(x, pred$mean, col="blue")
  lines(x, pred$mean - 2*sqrt(pred$sd2), col="blue") # nugs???
  lines(x, pred$mean + 2*sqrt(pred$sd2), col="blue") # nugs???
  legend(x = "topright", legend = c("data","truth", "wt avg", "95% UQ"),
         col = c("gray","black","red","blue"), lty = c(1,1,2,1))
}
hetgp_mse <- mean((pred$mean - y_true)^2)
hetgp_logs <- sum(logs.numeric(drop(y_true), family = "normal", mean = pred$mean, 
                             sd = sqrt(pred$sd2)))
toc <- proc.time()[3]
time3 <- toc - tic

# Store results ---------------------------------------------------------------

filename <- paste0("results/sims_", func, "_", setting, "_", r, ".csv")
if (file.exists(filename)) {
  results <- read.csv(filename)
  results <- rbind(results, c(seed, dgp_mse, dgp_r_mse,
                              deepgp_mse, hetgp_mse,
                              time1, time_r, time2, time3))
} else {
  results <- data.frame(seed = seed, dgp = dgp_mse, dgp_r = dgp_r_mse,
                        deepgp = deepgp_mse, hetgp = hetgp_mse,
                        time1 = time1, time_r = time_r,
                        time2 = time2, time3 = time3)
}
write.csv(results, filename, row.names = FALSE)

