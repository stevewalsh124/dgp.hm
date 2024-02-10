
library(dgp.hm)

library(MASS) #ginv
library(fields) #image.plot
library(mvtnorm) #rmvnorm
library(plgp) #distance (which is squared distances)
library(zoo) #rollmean
library(Matrix)

# Fit a shallow (one_layer) GP, or a (two_layer) deep GP (DGP)?
one_layer <- T

# Model the correlated errors with a covariance function?
# i.e., model Sigma_epsilon with a Matern/Exp2 cov fn?
err_cov <- "matern"#"exp2"#
loess_span <- 0.15 #LOESS-smoothing param for the mean
err_v   <- paste0("est", loess_span)#ifelse(err_cov == "matern", 2.5, 999)
err_g   <- NULL#sqrt(.Machine$double.eps)#
err_g_msg <- ifelse(is.null(err_g),"estg","fixg")

# Model the true power matter spectrum with cov fn (matern vs exp2)
cov_fn <- "matern"#"exp2"#

# Number of low res runs to avg over
nrun <- 16

# MCMC settings
nmcmc <- 1500#0
nburn <- 1000#0
kth <- 4

# Load data for a particular model --------------------------------------------

# Model: choose from 1-111 for training set, or c(0, 112:116) for test set
mte <- 1
if (mte <= 111) {
  mte_name <- paste0("M", if(mte < 100) {"0"}, if (mte < 10) {"0"}, mte) 
} else {
  test_names <- c("E001", "E002", "E003", "E009", "E010")
  mte_name <- test_names[mte - 111]
}

# first column is k, 2nd is linear pert theory, 3:18 is low-res, 19 is hi-res
file_name <- paste0("../Mira-Titan-IV-data/Mira-Titan-2021/STEP499/pk_", 
                    mte_name, "_test.dat")
pk2 <- read.table(file_name)
n <- nrow(pk2)

# Precisions/Inputs -----------------------------------------------------------

# Load precision info (k, prec_highres, prec_lowres, index_list)
load("../Mira-Titan-IV-Data/precision_and_indexes.Rdata")

# log10 of wavenumber (k) is X
x <- log10(k)

# Get precision info for eadh data type (low-res, hi-res, pert theory)
precs_lo <- ifelse(1:n %in% index_list$lowres.ix, prec_lowres, 0) * nrun
precs_hi <- ifelse(1:n %in% index_list$highres.ix, prec_highres, 0)
precs_pt <- ifelse(1:n %in% index_list$pert.ix, 10000, 0)

# Smooth the precision information (to remove steps)
Lam_lo <- rollmean(precs_lo, k = 10, fill = "extend")
Lam_hi <- rollmean(precs_hi, k = 10, fill = "extend")
Lam_pt <- rollmean(precs_pt, k = 10, fill = "extend")

# Precision matrix for weighted average, mu_z
Lam_z <- Lam_pt + Lam_lo + Lam_hi
precs <- diag(Lam_z)

# Response --------------------------------------------------------------------

# Grab responses
y_pt <- pk2[, 2]
y_low <- as.matrix(pk2[, 3:18])
y_hi <- pk2[, 19]

# Scale all of these to the emulation space (scrP function)
y_pt <- scrP(y_pt, k)
y_low <- scrP(y_low, k)
y_hi <- scrP(y_hi, k)

# Obtain low res average
y_lra <- rowMeans(y_low)

# get a weighted average (mu_z) across low, high and pert theory
# See "Weighted average from multiple computer experiments"
# in Walsh dissertation: 3.7.1 Appendix E: Derivations
mu_z <- (1 / Lam_z) * (Lam_pt * y_pt + Lam_lo * y_lra + Lam_hi * y_hi)

# STOPPED HERE FOR NOW 

# hi res precs on the low res index (for comparison)
precs_hi <- prec_highres
if(use_hi){
  # Weighted average of low res runs and hi-res
  hi_wt <- unique(prec_highres/prec_lowres)[1]
  y_avg <- mu_z
  # 16 low-res, and hi-res is 3.72x more precise than low-res
  # adjust effective sample size accordingly
  nrunn <- nrun + hi_wt
} else {
  y_avg <- y_lra
  nrunn <- nrun
}

# standardize the inputs (x) and outputs (y)
# get x to live in [0,1]
x <- (x - min(x))/(max(x)-min(x))

# make each y approx. mean 0, sd 1
mean_sz <- mean(y_avg)
sd_sz <- sd(y_avg)
y_avg <- (y_avg - mean_sz)/sd_sz
y_lra <- (y_lra - mean_sz)/sd_sz
y_hi <- (y_hi - mean_sz)/sd_sz
for (i in 1:nrow(Y)) { Y[i,] <- (Y[i,] - mean_sz)/sd_sz  }

# get (squared) distance matrix for the inputs
D <- plgp:::distance(x)

# adjust the precision info based on this standardization (sd_sz above)
precc <- precs*sd_sz^2
sdd <- sqrt(1/precc)
A <- diag(sdd)

# one way to get Sigma_hat
if(use_hi){
  n <- ncol(Y)
  r <- nrow(Y)
  
  sum_YYt <- matrix(0, n, n)
  for (i in 1:r) { sum_YYt <- sum_YYt + (Y[i,]-y_avg)%*%t(Y[i,]-y_avg) }
  sum_YYt <- sum_YYt + hi_wt*(y_hi-y_avg)%*%t(y_hi-y_avg)
  
  Sigma_hat <- 1/(r+hi_wt-1) * sum_YYt
} else {
  Sigma_hat <- cov(Y)
}

# This is the way to get Sigma_hat we will focus on
if(cf_errors){
  # # logl_cov* files use the same names (eg: logl_SW, fit_two_layer_SW)
  # if(!one_layer) source("logl_cov_1L.R")
  # varvec <- 1/precc
  # Sigma_hat_ho <- get_matern(x, Y_sim - colMeans(Y_sim), nmcmc = err_mcmc, nburn = err_burn, 
  #                         cov = err_cov, v = err_v, true_g = err_g, varvec = varvec)
  # Sigma_hat <- diag(sqrt(varvec)) %*% Sigma_hat_ho %*% diag(sqrt(varvec))
  # if(!one_layer) source("logl_cov.R")
  # D <- plgp::distance(x)
  if(loess_span != 0){
    loess_fit <- loess(y_avg ~ x, span = loess_span)
    avg_loess <- loess_fit$fitted
  } else {
    avg_loess <- y_avg
  }
  
  Y_zm <- matrix(NA, nrow(Y), ncol(Y))
  for (ii in 1:nrow(Y_zm)) Y_zm[ii,] <- Y[ii,] - avg_loess
  counter <- 0
  
  nl_matern <- function(par, D, Y, A) {
    theta <- exp(par[1])
    g <- exp(par[2])
    tau2 <- exp(par[3])
    kappa <- par[4]
    if(counter %% 50 == 0) print(c(theta,g,tau2,kappa))
    n <- ncol(Y)
    K <- tau2 * (geoR::matern(sqrt(D),phi = theta,kappa = kappa) + diag(g,n))
    
    Ki <- solve(K)
    ldetK <- determinant(K, logarithm=TRUE)$modulus
    
    ll <- 0
    for (ii in 1:nrow(Y)) {
      yit <- solve(A) %*% Y[ii,]
      ll <- ll - (1/2)*(t(yit) %*% Ki %*% yit) - (1/2)*ldetK 
    }
    counter <<- counter + 1
    Kii <<- Ki
    return(-ll)
  }
  
  lo_ind <- index_list$lowres.ix
  hi_ind <- index_list$highres.ix
  pt_ind <- index_list$pert.ix
  hi_only <- hi_ind[which(!(hi_ind %in% lo_ind))]
  
  out <- optim(c(.1, -5, .1, 1), nl_matern, method="L-BFGS-B", lower=c(-5,-35,-5,0.2),
               upper=c(4,3,4,100), D=D[lo_ind, lo_ind], Y=Y_zm[,lo_ind], A=A[lo_ind,lo_ind])
  
  # out <- optim(c(log(theta_e_true), log(g_e_true), log(tau2_e_true), true_v), 
  #              nl_matern, method="L-BFGS-B", lower=c(-5,-15,-5,0.2),
  #              upper=c(5,5,5,40), D=D, Y=Y_sim, A=A)
  
  kappa_hat <- out$par[4]
  phi_hat <- exp(out$par[1])
  theta_hat <- phi_hat*sqrt(kappa_hat)
  g_hat <- exp(out$par[2])
  tau2_hat <- exp(out$par[3])
  
  print(rbind(c("theta_e", "g_e", "tau2_e", "smooth"),
              c(theta_hat, g_hat, tau2_hat, kappa_hat)))
  
  Matern_hat <- tau2_hat * (geoR::matern(sqrt(D[lo_ind,lo_ind]), phi = phi_hat, kappa = kappa_hat) + diag(g_hat,ncol(Y[,lo_ind])))
  Sigma_hatty <-  A[lo_ind,lo_ind] %*% Matern_hat %*% A[lo_ind,lo_ind]
  # doubt this is right... need to redo the A matrices and precision info vectors
  Sigma_hat <- as.matrix(bdiag(solve(A[pt_ind,pt_ind])^0.5 %*% diag(rep(1/10000, length(pt_ind))) %*% solve(A[pt_ind,pt_ind])^0.5,
                               Sigma_hatty,
                               solve(A[hi_only,hi_only])^0.5 %*% diag((precs_hi[hi_only]*sd_sz^2)^-1) %*% solve(A[hi_only,hi_only])^0.5)) 
}

# if tapering Sigma_hat, use the bohman correlation function
if(taper_cov) Sigma_hat <- Sigma_hat * bohman(sqrt(D), tau=tau_b)

# plot the weighted average and also look at the values in Sigma_hat
par(mfrow=c(1,2))
matplot(x,t(Y),type="l",col="gray", ylim = range(y_hi, Y, y_avg, y_pt))
lines(x,y_hi,col="blue",lwd=2)
lines(x,y_avg,col="red",lwd=2)
legend("bottomright", legend = c("hi","avg"), lty=1, col=c("blue","red"))
image.plot(Sigma_hat/nrunn, main = "input as sigma_hat")

####################
# run for sim data #
####################

if(one_layer){
  fitcov <- fit_one_layer_SW(x = x, y = c(y_avg), nmcmc = nmcmc, true_g = 1e-10,
                             Sigma_hat = Sigma_hat/nrunn, cov = cov_fn)
} else {
  if(force_id_warp){
    fitcov <- fit_two_layer_SW(x = x, y = c(y_avg), nmcmc = nmcmc, Sigma_hat = Sigma_hat/nrunn, cov = cov_fn, pmx = pmx, true_g = 1e-10,
                               vecchia = vecchia, settings = list(alpha =list(theta_w=1000), beta=list(theta_w=.0001/1000)))
  } else {
    fitcov <- fit_two_layer_SW(x = x, y = c(y_avg), nmcmc = nmcmc, w_0 = w0_from_mte1_50k, true_g = 1e-10,
                               Sigma_hat = Sigma_hat/nrunn, cov = cov_fn, pmx = pmx,
                               vecchia = vecchia)
  }
}



fitcov <- trim(fitcov, nburn, kth)
plot(fitcov)

v <- fitcov$v

par(mfrow=c(1,1))
fitcov <- est.true(fitcov)
plot.true(fitcov)
if(mte %in% 1:111){
  cosmicEmu <- read.csv(paste0("R/CosmicEmu/2022-Mira-Titan-IV/P_tot/orig_111/EMU",
                               mte-1,".txt"),sep="", header = F)
} else {
  if(mte==0) cosmicEmu <- read.csv("R/CosmicEmu/2022-Mira-Titan-IV/P_tot/test_6/EMU5.txt",
                                   sep="", header = F)
  if(mte==112) cosmicEmu <- read.csv("R/CosmicEmu/2022-Mira-Titan-IV/P_tot/test_6/EMU0.txt",
                                     sep="", header = F)
  if(mte==113) cosmicEmu <- read.csv("R/CosmicEmu/2022-Mira-Titan-IV/P_tot/test_6/EMU1.txt",
                                     sep="", header = F)
  if(mte==114) cosmicEmu <- read.csv("R/CosmicEmu/2022-Mira-Titan-IV/P_tot/test_6/EMU2.txt",
                                     sep="", header = F)
  if(mte==115) cosmicEmu <- read.csv("R/CosmicEmu/2022-Mira-Titan-IV/P_tot/test_6/EMU3.txt",
                                     sep="", header = F)
  if(mte==116) cosmicEmu <- read.csv("R/CosmicEmu/2022-Mira-Titan-IV/P_tot/test_6/EMU4.txt",
                                     sep="", header = F)
  
}

cosmscrP <- log10(cosmicEmu[,2]*k^1.5/(2*pi^2))
cosmscrPsz <- (cosmscrP - mean_sz)/sd_sz
lines(fitcov$x, cosmscrPsz, col="green", lwd=1)

# Model 106 is hell
cb_cols   <- c("#999999", "#009E73", "#56B4E9", "#D55E00",
               "#E69F00", "#F0E442", "#0072B2", "#CC79A7")

# plot log10 space
# png(paste0("png/model",mte,"_emuspace.png"), width=4000, height = 2400, res=400)
par(mar=c(4,4.5,1,1), mfrow=c(1,1))
plot(fitcov$x, fitcov$y, type="n", #xlim = log10(c(.04,.35)),
     ylim = range(Y),
     xlab=expression(log[10](k)),
     ylab='script P')#TeX(r'($log_{10}(k^{1.5}P(k)/2\pi^2)$)'))
for (i in 1:16) lines(fitcov$x, Y[i,], col="gray", lwd=1)
lines(fitcov$x, colMeans(Y), col=cb_cols[2], lwd=2)
lines(fitcov$x, y_hi, col=cb_cols[3], lwd=2)
# lines(fitcov$x, fitcov$Ms[,1], col="red", lwd=2)
lines(fitcov$x, fitcov$ubb, col=cb_cols[4], lwd=2, lty=3)
lines(fitcov$x, fitcov$lbb, col=cb_cols[4], lwd=2, lty=3)
lines(fitcov$x, cosmscrPsz, col=cb_cols[8], lwd=2, lty=2)
legend("bottomright",legend = c("low res","low res avg","hi res", "UQ","cosmicEMU"), 
       col = c("gray",cb_cols[c(2,3,4,8)]), lty=c(1,1,1,3,2), lwd=2)
# dev.off()

m <- fitcov$m

# png(paste0("png/model",mte,"_emuspace_rmAvg.png"), width=4000, height = 2400, res=400)
par(mar=c(4,4.5,1,1), mfrow=c(1,1))
plot(log10(k), fitcov$y - m, type="n", #xlim = log10(c(.04,.35)),
     ylim = c(-.33,.33),
     xlab=expression(log[10](k)),
     ylab='script P')#TeX(r'($log_{10}(k^{1.5}P(k)/2\pi^2)$)'))
for (i in 1:16) lines(log10(k), Y[i,] - m, col="gray", lwd=1)
lines(log10(k), colMeans(Y) - m, col=cb_cols[2], lwd=2)
lines(log10(k), y_hi - m, col=cb_cols[3], lwd=2)
# lines(log10(k), fitcov$Ms[,1], col="red", lwd=2)
lines(log10(k), fitcov$ubb - m, col=cb_cols[4], lwd=2, lty=3)
lines(log10(k), fitcov$lbb - m, col=cb_cols[4], lwd=2, lty=3)
lines(log10(k), cosmscrPsz - m, col=cb_cols[8], lwd=2, lty=2)
legend("topleft",legend = c("low res","low res avg","hi res", "UQ","cosmicEMU"), 
       col = c("gray",cb_cols[c(2,3,4,8)]), lty=c(1,1,1,3,2), lwd=2)
# dev.off()


if(!one_layer) plot.warp(fitcov)

if(PDF) dev.off()

toc <- proc.time()[3]

(timing <- toc - tic)

fitcov$Cs <- fitcov$Ms <- fitcov$Ss <- 0
rm(Sigma_hat); rm(Sigma_hatty); rm(sum_YYt)
rm(A); rm(D); rm(Lam_hi); rm(Lam_lo)
rm(Lam_pt); rm(Lam_z); rm(Lam_zi)

if(saveImage) save.image(file = paste0("rda/emuspace_",nmcmc,"_",one_layer,if(pmx){"_pmx"},
                                       if(cf_errors){paste0("_cfe",err_v,err_g_msg)},if(taper_cov){paste0("tpr",tau_b)},
                                       if(force_id_warp){"_fiw"},if(vecchia){"_vec"},"model",mte,"_",k_sm,".rda"))


fitcov <- fit_one_layer_SW(x, y_avg, nmcmc = 5000, true_g = 1e-6, Sigma_hat = Sigma_hat)
  # could fix true_g = 1e-6
plot(fitcov) # from deepgp package


