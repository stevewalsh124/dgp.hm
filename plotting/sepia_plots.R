###############################################################################
# This script visualizes and interpolates the CAMB model runs, and also 
# obtains principal components and their corresponding weights 
# These will then get used to predict spectra for holdout cases
#
# Models 1-32 are used here.
#
###############################################################################

JPG = F

tic <- proc.time()[3]
library(dgp.hm)
library(GPfit)
xs <- read.csv("../fitting/results/CAMB/xcamb_int.csv")$x

cbcols <- palette.colors(palette = "Okabe-Ito")

# All the Runs and Models are in units Mpc/h and h/Mpc.
# To be able to plot everything in one plot, we have to convert to Mpc or 1/Mpc.
# To convert k (column 1), we need to multiply it by hubble (to get k in units of 1/Mpc) 
# while the power spectrum (column 2) needs to be divided by h^3 (to get it in Mpc^3).

# load the training design (ignore 7th column, indicating des=1 or des=2)
X = read.csv("../CosmicEmu_etc/pow64/cambDesigns_32x6x2.csv")

# standardize designs to be in the unit hypercube [0,1]^6
# ignore the last column of X, which denotes des=1 or des=2
for (i in 1:(ncol(X)-1)) X[,i] = (X[,i] - min(X[,i]))/(max(X[,i])-min(X[,i]))

# train on design 1, test on design 2
#i_train = which(X$des == 1)
i_train = 1:64
nruns = length(i_train)
des_train = X[i_train, -7]

p = ncol(des_train)  # number of input parameters
# hubble values to adjust k and P(k)
h <- read.csv("../CosmicEmu_etc/pow64/cambDesigns_32x6x2.csv")$h

# read in the linear output
run <- i_train[1]

par(mfrow=c(1,1))
sim <- read.table(paste0("../CosmicEmu_etc/pow64/RUN",run,"/oneh_matterpower.dat"))
k_c <- sim$V1 * h[run]
pk_c <- sim$V2/ (h[run])^3

# create k values of CAMB for interpolation
k_ci <- 10^xs

# interpolate to common k grid (k_ci)
pk_ci <- approx(k_c, scrP(pk_c, k_c), k_ci)$y

# Collect training runs
# So, dim of pk_cis is (# of interp k vals) x (# of CAMB runs)
pk_cs <- matrix(NA, length(pk_c), nruns)
pk_cis <- matrix(NA, length(pk_ci), nruns)
k_cs <- matrix(NA, length(k_c), nruns)
pk_cis[,1] <- pk_ci
pk_cs[,1] <- pk_c
k_cs[,1] <- k_c

for (i in 2:nruns) {
  run = i_train[i]
  sim <- read.table(paste0("../CosmicEmu_etc/pow64/RUN",
                           run,"/oneh_matterpower.dat"))
  k_cs[,i] <- k_c <- sim$V1 * h[run]
  pk_cs[,i] <- pk_c <- sim$V2 / (h[run])^3
  # lines(log10(k_c), scrP(pk_c, k_c), type="l", col=run)
  pk_cis[,i] <- approx(k_c, scrP(pk_c, k_c), k_ci)$y
  
}



preds_dgp <- read.table("../pca/dgp_camb_preds.txt")
preds_inf <- read.table("../pca/inf_res_preds.txt")


means <- apply(pk_cis,1,mean)
sds <- apply(pk_cis,1,sd)


mse.fun <- function(x,y) mean((x-y)^2)
mse.dgp <- mse.inf <- numeric(32)
for(i in 1:32){
  mse.dgp[i] <- mse.fun(pk_cis[,32+i], as.numeric(preds_dgp[i,]))
  mse.inf[i] <- mse.fun(pk_cis[,32+i], as.numeric(preds_inf[i,]))
}

mm <- apply(pk_cis[,1:32], 1, mean)
ss <- apply(pk_cis[,1:32], 1, sd)
llim <- -3.5; ulim <- 3.5

###################################################################################
# Predictive plot: DGP.FCO-trained vs inf-res-trained, compared to actual inf-res #
###################################################################################

if(JPG) jpeg("../paper/pred_diffs_CAMB.jpeg", width = 12, height = 6, units = "in", res=300)

par(mfrow=c(4,8),
    mar=c(0,0,0,0),
    oma=c(4,4,0.5,0))

yvals_all <- c()
for(i in 1:32){
  yvals_all <- c(yvals_all,
                 (preds_dgp[i,]-mm)/ss - (pk_cis[,32+i]-mm)/ss,
                 (preds_inf[i,]-mm)/ss - (pk_cis[,32+i]-mm)/ss)
}
ylims <- range(yvals_all, na.rm=TRUE)

for(i in 1:32){
  if(i > 25){
    plot(log10(k_ci), (pk_cis[,32+i]-mm)/ss-(pk_cis[,32+i]-mm)/ss,type="l", xlab=paste0(i),
         ylab="Script P", lty=2, yaxt="n", ylim=ylims, cex.axis=1.4)
  } 
  if(i %in% c(1,9,17)){
    plot(log10(k_ci), (pk_cis[,32+i]-mm)/ss - (pk_cis[,32+i]-mm)/ss,type="l", xlab=paste0(i),
         ylab="Script P", lty=2,
         xaxt="n", ylim=ylims, cex.axis=1.4)
  }
  if(i == 25){
    plot(log10(k_ci), (pk_cis[,32+i]-mm)/ss - (pk_cis[,32+i]-mm)/ss,type="l", xlab=paste0(i),
         ylab="Script P", lty=2, ylim=ylims, cex.axis=1.4)
  }
  if(i <= 24 & !(i %in% c(1,9,17,25))){
    plot(log10(k_ci), (pk_cis[,32+i]-mm)/ss - (pk_cis[,32+i]-mm)/ss,type="l", xlab=paste0(i),
         ylab="Script P", lty=2,
         xaxt="n",yaxt="n", ylim=ylims, cex.axis=1.4)
  }
  
  lines(xs,(preds_dgp[i,]-mm)/ss - (pk_cis[,32+i]-mm)/ss,col=cbcols[3], lwd=2, lty=1)
  lines(xs, (preds_inf[i,]-mm)/ss - (pk_cis[,32+i]-mm)/ss, col=cbcols[7], lwd=2, lty=3)
}
mtext(expression(log[10](k)), side = 1, outer = TRUE, line = 3)
mtext(expression(paste("\U1D4AB","  (k)", ", centered and scaled")), side = 2, outer = TRUE, line = 2)

if(JPG) dev.off()

####################################################
# Dot plot comparing MSE or RMSE of dgp vs inf_res #
####################################################

if(JPG) jpeg("../paper/mse_dot.jpeg", width = 3, height = 3, units = "in", res=300)
par(mfrow=c(1,1), mar=c(2,2,0,0),oma=c(2,2,0.5,0.5), pty="s")
plot(mse.dgp, mse.inf, col = cbcols[4], xlim = range(mse.dgp, mse.inf), 
     ylim = range(mse.dgp, mse.inf), xlab="", ylab="")
abline(0,1)
mtext(expression(MSE[DGP.FCO]), side = 1, outer = TRUE, line = 0.5)
mtext(expression(MSE[INF]), side = 2, outer = TRUE, line = 0.5)
if(JPG) dev.off()
