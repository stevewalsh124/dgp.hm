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
lines(x, f2(x, m2=m2,u2=u2,sd = 0))
lines(x, y_avg2, col="red", lty=2)
lines(x, fit$m, col="blue")
lines(x, fit$ub, col="blue")
lines(x, fit$lb, col="blue")
legend(x = "topright", legend = c("data","truth", "wt avg", "95% UQ"),
       col = c("gray","black","red","blue"), lty = c(1,1,2,1))

# repeat for r replicates (don't plot these)
mses1 <- c()
for (j in 1:r) {
  print(j)
  m1 <- runif(1, 0.5, 1.5)
  u1 <- runif(1, 1.5, 2.5)
  y1s <- matrix(nrow = r, ncol = n)
  for (i in 1:r) y1s[i,] <- f1(x, m1=m1, u1=u1)
  y_avg1 <- colMeans(y1s)
  var_y1 <- mean(apply(y1s, 2, var))
  Sigma1 = diag(var_y1, n)
  fit <- dgp.hm::fit_two_layer_hm(x, y_avg1, Sigma_hat = Sigma1, nmcmc = 7500)
  fit <- trim(fit, 2500, 5)
  fit <- est_true(fit)
  mses1[j] <- mean((fit$m - f1(x, m1=m1,u1=u1,sd = 0))^2)
}

mean(mses1)

# Boxplot
boxplot(mses1, ylim=c(0,.018), main="mses1")
grid(nx = NULL, ny = NULL,
     col = "#ebebeb", lwd = 2, lty=1)
boxplot(mses1, add = TRUE)

mses2 <- c()
for (j in 1:r) {
  print(j)
  m2 <- runif(1, 0.6, 1.4)
  u2 <- runif(1, 0.6, 1.4)
  y2s <- matrix(nrow = r, ncol = n)
  for (i in 1:r) y2s[i,] <- f2(x, m2=m2, u2=u2)
  y_avg2 <- colMeans(y2s)
  var_y2 <- mean(apply(y2s, 2, var))
  Sigma2 = diag(var_y2, n)
  fit <- dgp.hm::fit_two_layer_hm(x, y_avg2, Sigma_hat = Sigma2, nmcmc = 7500)
  fit <- trim(fit, 2500, 5)
  fit <- est_true(fit)
  mses2[j] <- mean((fit$m - f2(x, m2=m2,u2=u2,sd = 0))^2)
}

mean(mses2)

# Boxplot
boxplot(mses2, ylim=c(0,.005), main="mses2")
grid(nx = NULL, ny = NULL,
     col = "#ebebeb", lwd = 2, lty=1)
boxplot(mses2, add = TRUE)
