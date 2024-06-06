###############################################################################
# This script fits the Bayesian hierarchical model for simulated data under
# two different functions: f1 and f2
###############################################################################

library(dgp.hm)
library(mvtnorm)

# Simulated toy example from DPC paper
x <- seq(0, 4, by=0.1)
n <- length(x)

model <- 1 # set to 1 for "function 1", 2 for "function 2"
setting <- 4 #set to 1 for "Setting A", 2 for "B", 3 for "C"; 4 for "D"
if(!(model %in% 1:2)) stop("model should be 1 or 2")
if(!(setting %in% 1:4)) stop("setting should be 1, 2, 3, or 4")

r <- 5       # number of function realizations
n_sims <- 50 # number of sim studies to replicate

use_true_Sigma <- T # if T, use Sigma_true; if F, estimate Sigma_hat

args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))

# parameters for f1 and f2
m1 <- runif(1, 0.5, 1.5)
m2 <- runif(1, 0.6, 1.4)
u1 <- runif(1, 1.5, 2.5)
u2 <- runif(1, 0.6, 1.4)

# specify variance structure for the function
if(model == 1){
  if(setting == 1){
    sdvec <- rep(0.1, n)
    Sigma_true <- diag(sdvec^2)
  } 
  if(setting == 2){
    sdvec <- 1/(10*(x+1))
    Sigma_true <- diag(sdvec^2)
  } 
  if(setting == 3){
    sdvec <- 1/(10*(x+1))
    sdvec <- ifelse(x<=0.4, 5e-4, sdvec)
    sdvec <- ifelse(x>0.4 & x<3, 0.75*sdvec, sdvec)
    sdvec <- ifelse(x>=3, 0.5*sdvec, sdvec)
    Sigma_true <- diag(sdvec^2)
  }
  if(setting == 4){
    D <- deepgp::sq_dist(x)
    Sigma_true = deepgp:::Matern(distmat = D, tau2 = 0.1^2, 
                                 theta = 1e-2, g = 1e-8, v = 2.5)
  }
}

if(model == 2){
  if(setting == 1){
    sdvec <- rep(0.15, n)
    Sigma_true <- diag(sdvec^2)
  } 
  if(setting == 2){
    sdvec <- sqrt(abs(x-2))/10
    Sigma_true <- diag(sdvec^2)
  } 
  if(setting == 3){
    sdvec <- sqrt(abs(x-2))/10
    sdvec <- ifelse(x<=0.5, 5e-4, sdvec)
    sdvec <- ifelse(x>0.5 & x<2, 0.75*sdvec, sdvec)
    sdvec <- ifelse(x>=2, 0.5*sdvec, sdvec)
    Sigma_true <- diag(sdvec^2)
  }
  if(setting == 4){
    D <- deepgp::sq_dist(x)
    Sigma_true = deepgp:::Matern(distmat = D, tau2 = 0.15^2, 
                                 theta = 1e-2, g = 1e-8, v = 2.5)
  }
}

# make sure the sd plot looks right (compared to appropriate values)
image(Sigma_true)

# functions f1 and f2 to generate sims
f1 <- function(x, m1, u1, Sigma=diag(0,n,n)){
  w <- sqrt(25 - (u1/2)^2)
  return(m1*exp(-u1*x/2)*cos(w*x)-m1*x/5 + rmvnorm(1, sigma = Sigma))
}

f2 <- function(x, m2, u2, Sigma=diag(0,n,n)){
  return(exp(-m2*(x-3)^2) + exp(-u2*(x-1)^2) - 0.05*sin(8*(x-1.9)) + 
           rmvnorm(1, sigma = Sigma))
}

# Get n function evaluations for f1 and f2 
Y <- matrix(nrow = r, ncol = n)
if(model == 1){
  for (i in 1:r) Y[i,] <- f1(x, m1=m1, u1=u1, Sigma = Sigma_true)
} else {
  for (i in 1:r) Y[i,] <- f2(x, m2=m2, u2=u2, Sigma = Sigma_true)
}

# Get average and sigma
y_avg <- colMeans(Y)
if(model == 1) y_true <- f1(x, m1=m1, u1=u1)
if(model == 2) y_true <- f2(x, m2=m2, u2=u2)

# Plot these realizations
matplot(x, t(Y), ylab="f(x)", type="l")

# use either true cov structure, or estimate it
if(use_true_Sigma){
  Sigma_hat <- Sigma_true
} else {
  if(model %in% 1:3){
    var_y <- mean(apply(Y, 2, var))
    Sigma_hat = diag(var_y, n)
  } else {
    params_hat <- opt_matern(D, t(Y), sdd =rep(1,n))
    Sigma_hat <- deepgp:::Matern(D, params_hat$tau2_hat, params_hat$theta_hat,
                                 g=1e-8, v=2.5)
  }
}

# Fit DGP model, with true Sigma
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
  if(model == 1) y_true <- f1(x, m1=m1, u1=u1)
  if(model == 2) y_true <- f2(x, m2=m2, u2=u2)

  # generate draws from the true function, and get average
  Y <- matrix(nrow = r, ncol = n)
  if(model==1) for (i in 1:r) Y[i,] <- f1(x, m1=m1, u1=u1, Sigma = Sigma_true)
  if(model==2) for (i in 1:r) Y[i,] <- f2(x, m2=m2, u2=u2, Sigma = Sigma_true)
  y_avg <- colMeans(Y)
  
  # use either true cov structure, or estimate it
  if(use_true_Sigma){
    Sigma_hat <- Sigma_true
  } else {
    if(model %in% 1:3){
      var_y <- mean(apply(Y, 2, var))
      Sigma_hat = diag(var_y, n)
    } else {
      params_hat <- opt_matern(D, t(Y), sdd =rep(1,n))
      Sigma_hat <- deepgp:::Matern(D, params_hat$tau2_hat, params_hat$theta_hat,
                                   g=1e-8, v=2.5)
    }
  }
  
  # fit Deep GP hierarchical model
  fit <- dgp.hm::fit_two_layer_hm(x, y_avg, Sigma_hat = Sigma_hat, nmcmc = 7500)
  fit <- trim(fit, 2500, 5)
  fit <- est_true(fit)
  mses[j] <- mean((fit$m - y_true)^2)
}

mean(mses)

# Boxplot
boxplot(mses, ylim=c(0,ifelse(model==1,.018,.005)), main="mses")
grid(nx = NULL, ny = NULL,
     col = "#ebebeb", lwd = 2, lty=1)
boxplot(mses, add = TRUE)

write.csv(mses, file = paste0("results/sims_",model,"_",setting,"_"
                              ,r,"_",n_sims,".csv"))