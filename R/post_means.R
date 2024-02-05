# get the posterior mean from each DGP fit
# DGP fits done from 1D_real_study_full_emuspace.R

tic <- proc.time()[3]

vecchia <- F
pmx <- F
one_layer <- F
force_id_warp <- F

# Smooth the precision info so it's not a step function for each data type (pt, lo, hi)?
smooth_precs <- T
if(smooth_precs) k_sm <- 10 #rolling mean uses k numbers

# load the precision data (k, prec_highres, prec_lowres, index_list)
load("Mira-Titan-IV-data/precision_and_indexes.Rdata")

library(MASS) #ginv
library(fields) #image.plot
library(mvtnorm) #rmvnorm
library(plgp) #distance (which is squared distances)
library(zoo) #rollmean
library(Matrix)
# library(latex2exp) #TeX in axis label

# Taper the covariance matrix before the MCMC fit?
taper_cov <- F

# Do a kriging step?
krig <- F

# Use hi res in Ybar calculation?
use_hi <- T

# Model the correlated errors with a covariance function?
cf_errors <- T
if(cf_errors){
  err_cov <- "matern"#"exp2"#
  loess_span <- 0.15
  err_v   <- paste0("est", loess_span)#ifelse(err_cov == "matern", 2.5, 999)
  err_g   <- NULL#sqrt(.Machine$double.eps)#
  err_g_msg <- ifelse(is.null(err_g),"estg","fixg")
}

ncores <- 2
tolpower <- -10

cov_fn <- "matern"#"exp2"#

if(taper_cov) tau_b <- .2
nrun <- 16
nmcmc <- 1500#0
nburn <- 1000#0
kth <- 4

bte <- 3 # cols 3-18 are low res


args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))

# number of models to evaluate 
n_mte <- 34 #111

# collect and save the posterior means from the training set
post_means <- post_means_unstz <- matrix(NA, n_mte, 351)

for(mte in 1:n_mte){
  print(mte)
  load(file = paste0("rda/emuspace_",nmcmc,"_",one_layer,if(pmx){"_pmx"},
                     if(cf_errors){paste0("_cfe",err_v,err_g_msg)},if(taper_cov){paste0("tpr",tau_b)},
                     if(force_id_warp){"_fiw"},if(vecchia){"_vec"},"model",mte,"_",k_sm,".rda"))
  
  post_means[mte,] <- fitcov$m
  post_means_unstz[mte,] <- fitcov$m * sd_sz + mean_sz
  
  if(mte==1) plot(fitcov$x, post_means_unstz[mte,], ylim = c(-1.88,1.88), type="l", col=mte)
  if(mte!=1) lines(fitcov$x, post_means_unstz[mte,], ylim = c(-1.88,1.88), type="l", col=mte)
}

save(post_means, file = "rda/post_means.rda")
save(post_means_unstz, file = "rda/post_means_unstz.rda")

# collect and save the posterior means from the test set
post_means_test <- post_means_unstz_test <- matrix(NA, 6, 351)

for(mte in c(0,112:116)){
  print(mte)
  load(file = paste0("rda/emuspace_",nmcmc,"_",one_layer,if(pmx){"_pmx"},
                     if(cf_errors){paste0("_cfe",err_v,err_g_msg)},if(taper_cov){paste0("tpr",tau_b)},
                     if(force_id_warp){"_fiw"},if(vecchia){"_vec"},"model",mte,"_",k_sm,".rda"))

  if(mte==112) id <- 1
  if(mte==113) id <- 2
  if(mte==114) id <- 3
  if(mte==115) id <- 4
  if(mte==116) id <- 5
  if(mte==0) id <- 6
  post_means_test[id,] <- fitcov$m
  post_means_unstz_test[id,] <- fitcov$m * sd_sz + mean_sz
  
  if(mte==0) plot(fitcov$x, post_means_unstz_test[id,], ylim=c(-1.88, 1.88), type="l", col=id)
  if(mte!=0) lines(fitcov$x, post_means_unstz_test[id,], ylim=c(-1.88, 1.88), type="l", col=id)
}

rownames(post_means_test) <- rownames(post_means_unstz_test) <- c("E001","E002","E003","E009","E010","M000")
save(post_means_test, file = "rda/post_means_test.rda")
save(post_means_unstz_test, file = "rda/post_means_unstz_test.rda")

# # compare my DGP post means with Jared's proc. conv. means
# conv <- read.csv("csv/z_zero.txt",sep="",header = F)
# matplot(fitcov$x,t(conv - post_means_unstz),type="l",col="gray")
# abline(h=0,col="red",lty=2)
# 
# conv_test <- read.csv("csv/z_zero_test.txt",sep="",header = F)
# matplot(fitcov$x,t(conv_test - post_means_unstz_test),type="l",col="gray")
# abline(h=0,col="red",lty=2)
