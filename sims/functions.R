
# This file contains the functions needed to get f1, f2, and Sigma_true for
# the sim study from the DPC paper

# Parameters (these are randomly sampled once sourced)
m1 <- runif(1, 0.5, 1.5)
m2 <- runif(1, 0.6, 1.4)
u1 <- runif(1, 1.5, 2.5)
u2 <- runif(1, 0.6, 1.4)

f1 <- function(x, m1, u1, Sigma=diag(0,n,n)){
  w <- sqrt(25 - (u1/2)^2)
  return(m1*exp(-u1*x/2)*cos(w*x)-m1*x/5 + rmvnorm(1, sigma = Sigma))
}

f2 <- function(x, m2, u2, Sigma=diag(0,n,n)){
  return(exp(-m2*(x-3)^2) + exp(-u2*(x-1)^2) - 0.05*sin(8*(x-1.9)) + 
           rmvnorm(1, sigma = Sigma))
}

get_Sigma_true <- function(x, n, func = 1, setting = 1) {
  if(func == 1) {
    if(setting == 1) {
      sdvec <- rep(0.1, n)
      Sigma_true <- diag(sdvec^2)
    } else if(setting == 2) {
      sdvec <- 1/(10*(x+1))
      Sigma_true <- diag(sdvec^2)
    } else if(setting == 3) {
      sdvec <- 1/(10*(x+1))
      sdvec <- ifelse(x<=0.4, 5e-4, sdvec)
      sdvec <- ifelse(x>0.4 & x<3, 0.75*sdvec, sdvec)
      sdvec <- ifelse(x>=3, 0.5*sdvec, sdvec)
      Sigma_true <- diag(sdvec^2)
    } else if(setting == 4) {
      D <- deepgp::sq_dist(x)
      Sigma_true <- deepgp:::Matern(distmat = D, tau2 = 0.1^2, 
                                    theta = 1e-2, g = 1e-8, v = 2.5)
    } else if(setting == 5) {
      D <- deepgp::sq_dist(x)
      Sigma <- deepgp:::Matern(distmat = D, tau2 = 0.1^2, 
                               theta = 5e-2, g = 1e-8, v = 2.5)
      Sigma_true <- diag(sddtrue) %*% Sigma %*% diag(sddtrue)
    }
  } else if(func == 2) {
    if(setting == 1) {
      sdvec <- rep(0.15, n)
      Sigma_true <- diag(sdvec^2)
    } else if(setting == 2) {
      sdvec <- sqrt(abs(x-2))/10
      Sigma_true <- diag(sdvec^2)
    } else if(setting == 3) {
      sdvec <- sqrt(abs(x-2))/10
      sdvec <- ifelse(x<=0.5, 5e-4, sdvec)
      sdvec <- ifelse(x>0.5 & x<2, 0.75*sdvec, sdvec)
      sdvec <- ifelse(x>=2, 0.5*sdvec, sdvec)
      Sigma_true <- diag(sdvec^2)
    } else if(setting == 4) {
      D <- deepgp::sq_dist(x)
      Sigma_true <- deepgp:::Matern(distmat = D, tau2 = 0.15^2, 
                                   theta = 1e-2, g = 1e-8, v = 2.5)
    } else if(setting == 5) {
      D <- deepgp::sq_dist(x)
      Sigma <- deepgp:::Matern(distmat = D, tau2 = 0.1^2, 
                               theta = 5e-2, g = 1e-8, v = 2.5)
      Sigma_true <- diag(sddtrue) %*% Sigma %*% diag(sddtrue)
    }
  }
  return(Sigma_true)
}
