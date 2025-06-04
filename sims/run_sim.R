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

library(dgp.hm) # must install locally
library(dpc)
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
if(!(setting %in% 1:6)) stop("setting should be 1, 2, 3, 4, 5, or 6")

set.seed(seed)
source("functions.R") # generates random model parameters upon sourcing

r <- 5 # number of replicates observed
vis <- FALSE # should plots be generated

# Generate data ---------------------------------------------------------------

x <- seq(0, 4, by = 0.1)
n <- length(x)
if(setting==5) sddtrue <- 1/exp(seq(0,2,length=length(x)))
if(setting==6) sddtrue <- 1/1.5^(seq(0,2,length=length(x)))
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
} else if(setting==4) {
  # Get initial tau2 estimate
  tau2hat_0 <- mean(apply(Yzm, 2, var))
  
  estimate_theta <- function(x, max = TRUE) {
    # count the number of peaks (local maxima)
    if (max == FALSE) x <- x * (-1)
    res <- rep(FALSE, length(x))
    if (x[1] > x[2]) res[1] <- TRUE
    if (x[length(x)-1] < x[length(x)]) res[length(res)] <- TRUE
    for (i in (2:(length(x)-1))) {
      if ((x[i-1] < x[i]) & (x[i+1] < x[i])) res[i] <- TRUE
    }
    n_peaks <- sum(res)
    # from n_peaks, estimate theta
    theta_hat <- exp(3.649 - 2.695*log(n_peaks))
    if(log(n_peaks) < 1) theta_hat <- exp(1.22)
    if(log(n_peaks) > 5) theta_hat <- exp(-10)
    return(theta_hat)
  }
  
  # Get initial theta estimate
  thetahat_0 <- mean(apply(t(Yzm), 2, estimate_theta))
  
  dx <- sq_dist(x)
  # Optimize parameter estimates on transformed scale
  params_hat <- opt_matern(dx, t(Yzm), sdd = rep(1,n), 
                           init = c(log(tau2hat_0/thetahat_0), 
                                    log(tau2hat_0)))
  Sigma_hat <- deepgp:::Matern(dx, params_hat$tau2_hat, params_hat$theta_hat,
                               g = 1e-8, v = 2.5)
} else if(setting %in% 5:6) {
  # Get initial tau2 estimate
  tau2hat_0 <- mean(apply(Yzm, 2, var))
  
  estimate_theta <- function(x, max = TRUE) {
    # count the number of peaks (local maxima)
    if (max == FALSE) x <- x * (-1)
    res <- rep(FALSE, length(x))
    if (x[1] > x[2]) res[1] <- TRUE
    if (x[length(x)-1] < x[length(x)]) res[length(res)] <- TRUE
    for (i in (2:(length(x)-1))) {
      if ((x[i-1] < x[i]) & (x[i+1] < x[i])) res[i] <- TRUE
    }
    n_peaks <- sum(res)
    # from n_peaks, estimate theta
    theta_hat <- exp(3.649 - 2.695*log(n_peaks))
    if(log(n_peaks) < 1) theta_hat <- exp(1.22)
    if(log(n_peaks) > 5) theta_hat <- exp(-10)
    return(theta_hat)
  }
  
  # Get initial theta estimate
  thetahat_0 <- mean(apply(t(Yzm), 2, estimate_theta))
  
  dx <- sq_dist(x)
  # Optimize parameter estimates on transformed scale
  params_hat <- opt_matern(dx, t(Yzm), sdd = sddtrue, 
                           init = c(log(tau2hat_0/thetahat_0), 
                                    log(tau2hat_0)))
  Sigma_hat <- diag(sddtrue) %*% 
    deepgp:::Matern(dx, params_hat$tau2_hat, 
                    params_hat$theta_hat, g = 1e-8, v = 2.5) %*% diag(sddtrue)
}

fit <- dgp.hm::fit_two_layer_hm(x, y_avg, Sigma_hat = Sigma_hat, nmcmc = 15000)
fit <- dgp.hm::trim(fit, 5000, 1)
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
fit1a <- dgp.hm::fit_two_layer_hm(x, y_avg, Sigma_hat = Sigma_hat/r, 
                                  nmcmc = 15000)
if(vis) plot(fit1a)
fit1a <- dgp.hm::trim(fit1a, 5000, 1)
fit1a <- dgp.hm::est_true(fit1a, return_all = TRUE)

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
dgp_r_logs <- sum(logs.numeric(drop(y_true), family = "normal", mean = fit1a$m, 
                               sd = apply(fit1a$Ss, 1, sd)))
toc <- proc.time()[3]
dgp_r_time <- toc - tic

# deepgp model ----------------------------------------------------------------

tic <- proc.time()[3]
fit2 <- deepgp::fit_two_layer(x_all, y_all, nmcmc = 15000, vecchia = T)
if(vis) plot(fit2)
fit2 <- deepgp::trim(fit2, 5000, 1)
fit2 <- predict(fit2, x)
if(vis) plot(fit2)
deepgp_mse <- mean((fit2$mean - y_true)^2)
deepgp_logs <- sum(logs.numeric(drop(y_true), family = "normal", 
                                mean = fit2$mean, sd = sqrt(fit2$s2/n)))
toc <- proc.time()[3]
deepgp_time <- toc - tic

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
hetgp_logs <- sum(logs.numeric(drop(y_true), family = "normal", 
                               mean = pred$mean, sd = sqrt(pred$sd2)))
toc <- proc.time()[3]
hetgp_time <- toc - tic

# Deep Process Convolutions (DPC) results -------------------------------------
# For the DPC model, Y and vars should be a length nfuncs list.
# List element i holds the n_i sets of observations of function i.
# Here n_i = r.
tic <- proc.time()[3]
Y_list = list(); vars_list = list()
for(i in 1:1){ # if you only have one function
  Y_list[[i]] = list(); vars_list[[i]] = list()
  for(j in 1:r) {
    Y_list[[i]][[j]] = cbind(x, Y[j, ])
    vars_list[[i]][[j]] = diag(Sigma_true)
  }
}

# Run DPC model
n_u=30; n_v=10
dpc_out = dpc(Y=Y_list, vars=vars_list, n_u=n_u, n_v=n_v, 
              nmcmc=15000, burn=5000,
              sprop_tau2_u=0.005, sprop_delta=0.3,
              sprop_v=rep(0.1, n_v))

# Get the DPC function estimate, LL, and UL
out_f = plot.dpc(output=dpc_out, graph=FALSE, return_mean=TRUE, 
                 return_draws=TRUE, gridpoints=x)
f_es = out_f$pmean[,-1]
f_ll = apply(out_f$pdraws, c(1,3), function(x) quantile(x,0.025, na.rm=T))
f_ul = apply(out_f$pdraws, c(1,3), function(x) quantile(x,0.975, na.rm=T))
if(vis){
  plot(x_all, y_all)
  lines(x, f_es, col='blue')
  lines(x, f_ll, lty=3, col='blue')
  lines(x, f_ul, lty=3, col='blue')
}
dpc_mse <- mean((f_es - y_true)^2)
dpc_logs <- sum(logs.numeric(drop(y_true), family = "normal", 
                               mean = f_es, sd = sqrt(apply(out_f$pdraws, 1, var))))
toc <- proc.time()[3]
dpc_time <- toc - tic

# Store results ---------------------------------------------------------------

filename <- paste0("results/sims_", func, "_", setting, "_", r, ".csv")
if (file.exists(filename)) {
  results <- read.csv(filename)
  results <- rbind(results, 
                   c(seed, dgp_mse, dgp_r_mse, deepgp_mse, hetgp_mse, dpc_mse,
                     dgp_logs, dgp_r_logs, deepgp_logs, hetgp_logs, dpc_logs,
                     dgp_time, dgp_r_time, deepgp_time, hetgp_time, dpc_time))
} else {
  results <- data.frame(seed = seed, dgp_mse = dgp_mse, dgp_r_mse = dgp_r_mse,
                        deepgp_mse = deepgp_mse, hetgp_mse = hetgp_mse, 
                        dpc_mse = dpc_mse, dgp_logs = dgp_logs, 
                        dgp_r_logs = dgp_r_logs, deepgp_logs = deepgp_logs, 
                        hetgp_logs = hetgp_logs, dpc_logs = dpc_logs,
                        dgp_time = dgp_time, dgp_r_time = dgp_r_time,
                        deepgp_time = deepgp_time, hetgp_time = hetgp_time, 
                        dpc_time = dpc_time)
}

write.csv(results, filename, row.names = FALSE)

