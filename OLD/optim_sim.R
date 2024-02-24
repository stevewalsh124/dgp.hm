############################################
# Basic sim study to estimate tau2 & theta #
# true values of tau2/theta/sdd set to     #
# emulate the real data                    #
############################################

library(dgp.hm)

n = 351
r = 16
x <- seq(-3,1,length=n)

dx <- deepgp::sq_dist(x)

tau2 <- 50
theta <- .005
# approximate the precision info from Mira-Titan-IV
sdd <- 1/sqrt(exp(seq(8.03,16.77,length=n)))
covmtx <- diag(sdd) %*% deepgp:::Matern(dx, tau2, theta, 1e-8, 2.5) %*% diag(sdd)

nsims <- 10
hats <- matrix(NA, nsims, 2)

for (i in 1:nsims) {
  
  if(i%%5==0) print(i)
  
  Y <- matrix(NA,n,r)
  for (j in 1:r) Y[,j] <- y <- t(chol(covmtx)) %*% rnorm(n)
  
  # matplot(x,Y, type="l")
  params <- opt_matern(dx = dx, y = y, sdd = sdd, init=c(exp(2),exp(2)))
  
  hats[i,] <- c(thetahat=params$theta_hat, tau2hat=params$tau2_hat)
}

par(mfrow=c(1,2))
hist(hats[,1], main = "theta hats"); abline(v=theta, col="blue", lwd=2)
hist(hats[,2], main = "tau2 hats"); abline(v=tau2, col="blue", lwd=2)
var(hats[,1])
var(hats[,2])
