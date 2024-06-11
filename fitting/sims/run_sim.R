###############################################################################
# This script conducts a "bake-off" of various surrogate models under two
# different test functions (f1 and f2) for 50 Monte Carlo repetitions.
#
# Important Arguments:
#   func: 1 or 2
#   setting: 1, 2, 3, or 4 (indicates noise settings A-D from DPC paper)
#   r: number of function realizations
#
# Output:
#   data frame of MSE values for each model are written to "results" folder
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
source("functions.R")

# Get settings from command line
func <- 1
setting <- 1
r <- 5
args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))
if(!(func %in% 1:2)) stop("func should be 1 or 2")
if(!(setting %in% 1:4)) stop("setting should be 1, 2, 3, or 4")
cat("func is ", func, "\n")
cat("setting is ", setting, "\n")
cat("r is ", r, "\n")

vis <- FALSE # should plots be generated
n_sims <- 50 # number of sim studies to replicate
use_true_Sigma <- T # if T, use Sigma_true; if F, estimate Sigma_hat

# Generate data
x <- seq(0, 4, by=0.1)
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

# Get average and truth
y_avg <- colMeans(Y)
if(func == 1) y_true <- f1(x, m1=m1, u1=u1)
if(func == 2) y_true <- f2(x, m2=m2, u2=u2)

# Use either true cov structure, or estimate it (for dgp.hm model only)
if(use_true_Sigma) {
  Sigma_hat <- Sigma_true
} else {
  if(setting %in% 1:3) {
    var_y <- mean(apply(Y, 2, var))
    Sigma_hat = diag(var_y, n)
  } else {
    params_hat <- opt_matern(D, t(Y), sdd =rep(1,n))
    Sigma_hat <- deepgp:::Matern(D, params_hat$tau2_hat, params_hat$theta_hat,
                                 g=1e-8, v=2.5)
  }
}

# Model 1: dgp.hm
fit <- dgp.hm::fit_two_layer_hm(x, y_avg, Sigma_hat = Sigma_hat, nmcmc = 20000)
plot(fit)
fit <- trim(fit, 15000, 5)
plot(fit)
fit <- est_true(fit)

# plot estimated function alongside data and avg
matplot(x, t(Y), type="l", col="gray")
lines(x, y_true)
lines(x, y_avg, col="red", lty=2)
lines(x, fit$m, col="blue")
lines(x, fit$ub, col="blue")
lines(x, fit$lb, col="blue")
legend(x = "topright", legend = c("data","truth", "wt avg", "95% UQ"),
       col = c("gray","black","red","blue"), lty = c(1,1,2,1))

# repeat for r replicates (don't plot these)
mses <- c()
for (j in 1:n_sims) {
  # generate random model parameters
  print(j)
  m1 <- runif(1, 0.5, 1.5)
  u1 <- runif(1, 1.5, 2.5)

  # generate a true function
  if(func == 1) y_true <- f1(x, m1=m1, u1=u1)
  if(func == 2) y_true <- f2(x, m2=m2, u2=u2)

  # generate draws from the true function, and get average
  Y <- matrix(nrow = r, ncol = n)
  if(func==1) for (i in 1:r) Y[i,] <- f1(x, m1=m1, u1=u1, Sigma = Sigma_true)
  if(func==2) for (i in 1:r) Y[i,] <- f2(x, m2=m2, u2=u2, Sigma = Sigma_true)
  y_avg <- colMeans(Y)
  
  # use either true cov structure, or estimate it
  if(use_true_Sigma){
    Sigma_hat <- Sigma_true
  } else {
    if(func %in% 1:3){
      var_y <- mean(apply(Y, 2, var))
      Sigma_hat = diag(var_y, n)
    } else {
      params_hat <- opt_matern(D, t(Y), sdd =rep(1,n))
      Sigma_hat <- deepgp:::Matern(D, params_hat$tau2_hat, params_hat$theta_hat,
                                   g=1e-8, v=2.5)
    }
  }
  
  # fit Deep GP hierarchical func
  fit <- dgp.hm::fit_two_layer_hm(x, y_avg, Sigma_hat = Sigma_hat, nmcmc = 7500)
  fit <- trim(fit, 2500, 5)
  fit <- est_true(fit)
  mses[j] <- mean((fit$m - y_true)^2)
}

mean(mses)

# Boxplot
boxplot(mses, ylim=c(0,ifelse(func==1,.018,.005)), main="mses")
grid(nx = NULL, ny = NULL,
     col = "#ebebeb", lwd = 2, lty=1)
boxplot(mses, add = TRUE)

write.csv(mses, file = paste0("results/sims_",func,"_",setting,"_"
                              ,r,"_",n_sims,".csv"))