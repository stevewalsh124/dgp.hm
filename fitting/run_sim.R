###############################################################################
# This script fits the Bayesian hierarchical model for simulated data under
# two different functions: f1 and f2
###############################################################################

library(dgp.hm)

# Simulated toy example from DPC paper
x <- seq(0, 4, by=0.1)
n <- length(x)
r <- 50 # number of replicates

# parameters for f1 and f2
m1 <- runif(1, 0.5, 1.5)
m2 <- runif(1, 0.6, 1.4)
u1 <- runif(1, 1.5, 2.5)
u2 <- runif(1, 0.6, 1.4)

# functions f1 and f2 to generate sims
f1 <- function(x, m1, u1, sd=0.1){
  w <- sqrt(25 - (u1/2)^2)
  return(m1*exp(-u1*x/2)*cos(w*x)-m1*x/5 + rnorm(length(x), sd=sd))
}

f2 <- function(x, m2, u2, sd=0.15){
  return(exp(-m2*(x-3)^2) + exp(-u2*(x-1)^2) - 0.05*sin(x-1.9) + rnorm(length(x), sd=sd))
}

# Get n function evaluations for f1 and f2 
y1s <- y2s <- matrix(nrow = r, ncol = n)
for (i in 1:r) y1s[i,] <- f1(x, m1=m1, u1=u1)
for (i in 1:r) y2s[i,] <- f2(x, m2=m2, u2=u2)

# Get average and sigma
y_avg1 <- colMeans(y1s)
y_avg2 <- colMeans(y2s)
Sigma1 = diag(0.1^2, n)
Sigma2 = diag(0.15^2, n)

# Plot these realizations
matplot(x, t(y1s), ylab="f1(x)", type="l")
matplot(x, t(y2s), ylab="f2(x)", type="l")

# Fit DGP model to f1
fit <- dgp.hm::fit_two_layer_hm(x, y_avg1, Sigma_hat = Sigma1, nmcmc = 7500)
plot(fit)
fit <- trim(fit, 2500, 5)
plot(fit)
fit <- est_true(fit)

# plot estimated function alongside data and avg
xx <- seq(0,4,length=100)
matplot(x, t(y1s), type="l", col="gray")
lines(xx, f1(xx, m1=m1,u1=u1,sd = 0))
lines(x, y_avg1, col="red", lty=2)
lines(x, fit$m, col="blue")
lines(x, fit$ub, col="blue")
lines(x, fit$lb, col="blue")
legend(x = "topright", legend = c("data","truth", "wt avg", "95% UQ"),
       col = c("gray","black","red","blue"), lty = c(1,1,2,1))


# Fit DGP model to f2
fit <- dgp.hm::fit_two_layer_hm(x, y_avg2, Sigma_hat = Sigma2, nmcmc = 10000)
plot(fit)
fit <- trim(fit, 5000, 5)
plot(fit)
fit <- est_true(fit)

# plot estimated function alongside data and avg
xx <- seq(0,4,length=100)
matplot(x, t(y2s), type="l", col="gray")
lines(xx, f2(xx, m2=m2,u2=u2,sd = 0))
lines(x, y_avg2, col="red", lty=2)
lines(x, fit$m, col="blue")
lines(x, fit$ub, col="blue")
lines(x, fit$lb, col="blue")
legend(x = "topright", legend = c("data","truth", "wt avg", "95% UQ"),
       col = c("gray","black","red","blue"), lty = c(1,1,2,1))
