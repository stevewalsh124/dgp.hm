############################################
# Basic sim study to estimate tau2 & theta #
# true values of tau2/theta/sd_p set to     #
# emulate the real data                    #
############################################

library(dgp.hm)
library(fields) # image.plot

# n: number of observed points per model run
# r: number of model runs
n = 351
r = 16

# input x reflects log10(k) values
x <- read.csv("../OLD/x.csv")$x

# create matrix of squared distances for covariance
dx <- deepgp::sq_dist(x)

# set true param values near those estimated from Mira-Titan-IV
tau2 <- 50
theta <- .005

# use mean and  precision info based on Mira-Titan-IV
y_true <- read.csv("../OLD/y_true.csv")$x
prec <- read.csv("../OLD/prec.csv")$x

# standard deviation from precision info
sd_p <- 1/sqrt(prec)

# create nonstationary covariance using precision info
covmtx <- diag(sd_p) %*% deepgp:::Matern(dx, tau2, theta, 1e-8, 2.5) %*% diag(sd_p)

# how many simulation studies to run?
nsims <- 10

# keep track of estimated parameters
hats <- matrix(NA, nsims, 2)

for (i in 1:nsims) {
  
  if(i%%5==0) print(i)
  
  # create r observations
  Y <- matrix(NA,n,r)
  for (j in 1:r) Y[,j] <- y <- t(chol(covmtx)) %*% rnorm(n) + y_true
  
  # standardize: subtract LOESS mean and scale based on sd
  y_bar <- loess(rowMeans(Y)~x,span = .15)$fitted
  sd_y <- sd(y_bar)
  # matplot(x,Y, type="l")
  
  # adjust precision info (in form of sd_p) by the scaling (sd_y)
  sd_s <- sd_p/sd_y
  params <- opt_matern(dx = dx, y = (y - y_bar)/sd_y, sdd = sd_s, init=c(exp(2),exp(2)))
  
  # estimate hyperparameters and construct covariance matrix
  hats[i,] <- c(thetahat=params$theta_hat, tau2hat=params$tau2_hat)
  Sigma_hat <- diag(sd_s) %*% 
    deepgp:::Matern(dx,params$tau2_hat, params$theta_hat,1e-8,2.5) %*% diag(sd_s)
  
  # compare estimated and true (after scaling with sd_y) covariance matrices
  # par(mfrow=c(1,2))
  # fields::image.plot(Sigma_hat, main = "estimated covmtx: Sigma_hat", 
  #                    zlim=range(covmtx/(sd_y^2), Sigma_hat))
  # fields::image.plot(covmtx /(sd_y^2), main = "true covmtx", 
  #                    zlim=range(covmtx/(sd_y^2), Sigma_hat))
}

# compare parameter estimates with true values
par(mfrow=c(1,2))
hist(hats[,1], main = "theta hats"); abline(v=theta, col="blue", lwd=2)
hist(hats[,2], main = "tau2 hats"); abline(v=tau2, col="blue", lwd=2)
