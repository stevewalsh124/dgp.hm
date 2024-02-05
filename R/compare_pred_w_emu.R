# compare predictions... my DGP and GP-PC combo vs cosmic emu
PDF <- T

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


# Use the high res runs to train?
hi <- F

# Which design do you want to train with? a, c, or both?
atf <- T
ctf <- T

seed <- 1
des <- "test6"

n_pc <- 10
pwdExp <- 1.95

args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))

if(PDF) pdf(paste0("pdf/compare_pred_w_emu_",n_pc,"_",pwdExp,".pdf"))

mses <- mses_emu <- mses_conv <- matrix(NA, 6, 351)

# import Kelly and Jared's process convolution estimates
# convln_pred <- read.csv("csv/z_zero_test.txt", sep = "", header=F)

for(mtee in c(0,112:116)){
  print(mtee)
  load(file = paste0("rda/emuspace_",nmcmc,"_",one_layer,if(pmx){"_pmx"},
                     if(cf_errors){paste0("_cfe",err_v,err_g_msg)},if(taper_cov){paste0("tpr",tau_b)},
                     if(force_id_warp){"_fiw"},if(vecchia){"_vec"},"model",mtee,"_",k_sm,".rda"))
  
  # png(paste0("png/model",mte,"_emuspace_rmAvg.png"), width=4000, height = 2400, res=400)
  par(mar=c(4,4.5,1,1), mfrow=c(1,1))
  plot(log10(kl), fitcov$y - m, type="n", #xlim = log10(c(.04,.35)),
       ylim = c(-.12,.12),
       xlab=expression(log[10](k)),
       ylab='script P')#TeX(r'($log_{10}(k^{1.5}P(k)/2\pi^2)$)'))
  for (i in 1:16) lines(log10(kl)[index_list$lowres.ix], (Y[i,] - m)[index_list$lowres.ix], col="gray", lwd=1)
  # lines(log10(kl)[index_list$lowres.ix], (colMeans(Y) - m)[index_list$lowres.ix], col=cb_cols[2], lwd=2)
  lines(log10(kl)[index_list$highres.ix], (y_hi - m)[index_list$highres.ix], col=cb_cols[3], lwd=2)
  y_pt <- (y_pt - mean_sz)/sd_sz
  lines(log10(kl)[index_list$pert.ix], (y_pt - m)[index_list$pert.ix], col=cb_cols[2], lwd=3)
  # lines(log10(kl), fitcov$Ms[,1], col="red", lwd=2)
  lines(log10(kl), fitcov$ubb - m, col=cb_cols[4], lwd=2, lty=3)
  lines(log10(kl), fitcov$lbb - m, col=cb_cols[4], lwd=2, lty=3)
  lines(log10(kl), cosmscrPsz - m, col=cb_cols[8], lwd=2, lty=2)
  legend("topright",legend = c("low res","pert","hi res", "UQ","cosmicEMU","DGP-GPPC","conv"), 
         col = c("gray",cb_cols[c(2,3,4,8)],"black","red"), lty=c(1,1,1,3,2,2), lwd=2)
  mtext(mte)
  # dev.off()
  
  meanterm = m
  
  if(mtee==0) id <- 6
  if(mtee==112) id <- 1
  if(mtee==113) id <- 2
  if(mtee==114) id <- 3
  if(mtee==115) id <- 4
  if(mtee==116) id <- 5
  # convln_pred[id,] <- (convln_pred[id,] - mean_sz)/sd_sz
  
  # lines(log10(kl), convln_pred[id,] - m, col="red")
  
  load(paste0("rda/cosmo-SA-",
                    ifelse(hi,"hi-","lo-"),"train",if(atf){"a"},if(ctf){"c"},"-",
                    "pred",des,if(des%in%c("FF","FFS")){n_ff},
                    if(des%in%c("unif","lhs","lhsS","unifS")){m},
                    if(des %in% c("lhsS","unifS")){paste0("s",seed)},
                    "-nPC",n_pc,"_pwdExp",pwdExp,".rda"))
  
  lines(log10(kl), etaEmu[,id] - meanterm, col="black", lwd=1, lty=2)

  mses[id,] <- (etaEmu[,id] - meanterm)^2
  mses_emu[id,] <- (cosmscrPsz - meanterm)^2
  # mses_conv[id,] <- (as.numeric(convln_pred[id,]) - meanterm)^2
  cat(c("RMSE_dgp",sqrt(sum((etaEmu[,id] - meanterm)^2))),"\n")
  cat(c("RMSE_emu",(sqrt(sum((cosmscrPsz - meanterm)^2)))),"\n")
  # cat(c("RMSE_conv",sqrt(sum((convln_pred[id,] - meanterm)^2))),"\n")
  
}

plot(kvals,colMeans(sqrt(mses)), type="l",ylim = range(colMeans(sqrt(mses)),colMeans(sqrt(mses_emu))),
     main = "average RMSE by k value")
lines(kvals,colMeans(sqrt(mses_emu)), col=cb_cols[8])
lines(kvals,colMeans(sqrt(mses_conv)), col="red")
legend("topright",legend = c("cosmicEMU","DGP-GPPC","conv"), 
       col = c(cb_cols[8],"black","red"), lwd=2)

boxplot(c(sqrt(mses)),c(sqrt(mses_conv)),c(sqrt(mses_emu)), names=c("DGP-GPPC","conv","cosmicEMU"), ylab="RMSE",
        main = "all RMSEs across k values and cosmologies", col = c(0,cb_cols[8],"red"))

if(PDF) dev.off()
