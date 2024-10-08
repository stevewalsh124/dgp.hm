###############################################################################
# This script visualizes and interpolates the CAMB model runs, and also 
# obtains principal components and their corresponding weights 
# These will then get used to predict spectra for holdout cases
#
# Models 1-32 are used here.
#
###############################################################################

tic <- proc.time()[3]

library(GPfit)
library(dgp.hm)

plot.ratio <- T

PDF <- F
if(PDF) pdf(paste0("CAMB_pred_des1",if(plot.ratio){"ratio"},".pdf"))

i_test = which(X$des != 2)
des_test <- read.csv("../CosmicEmu_etc/pow64/cambDesigns_32x6x2.csv")[i_test,-7]
ntest <- nrow(des_test)

# All the Runs and Models are in units Mpc/h and h/Mpc.
# To be able to plot everything in one plot, we have to convert to Mpc or 1/Mpc.
# To convert k (column 1), we need to multiply it by hubble (to get k in units of 1/Mpc) 
# while the power spectrum (column 2) needs to be divided by h^3 (to get it in Mpc^3).

# contains hubble values to adjust k and P(k)
h <- read.csv("../CosmicEmu_etc/pow64/cambDesigns_32x6x2.csv")$h

# create k values of CAMB for interpolation
k_ci <- k_ci <- 10^(seq(-1.7, 0.7, by=0.015))

# Collect remaining runs 2-32
# So, dim of pk_cis is (# of interp k vals) x (# of CAMB runs)
pk_cs <- matrix(NA, 611, 32)
pk_cis <- matrix(NA, length(k_ci), 32)
k_cs <- matrix(NA, 611, 32)
MSEs <- MSEs_Delta <- rep(NA,32)

par(mfrow=c(1,1))
# Working with runs 1-32; hold out one of them
for(holdout in 1:ntest){
  run = i_test[holdout]

  # read in emulator
  etaEmu <- read.csv(paste0("../pca/etaEmu10_1.9CAMB","des1",".csv"))[,holdout]

  # read in heldout value
  sim <- read.table(paste0("../CosmicEmu_etc/pow64/RUN",
                           run,"/oneh_matterpower.dat"))
  
  # record k, P(k), P(k) interpolated, and MSEs
  k_cs[,holdout] <- k_c <- sim$V1 * h[run]
  pk_cs[,holdout] <- pk_c <- sim$V2 / (h[run])^3
  pk_cis[,holdout] <- approx(k_c, scrP(pk_c, k_c), k_ci)$y
  MSEs[holdout] <- mean((etaEmu - pk_cis[,holdout])^2)

  
  # Plot the interpolated CAMB runs, or a ratio of Delta
  if(plot.ratio){
    
    # convert the scriptP(k) to Delta^2(k)
    Delta_holdout <- sqrt((10^pk_cis[,holdout])*(k_ci^1.5))
    Delta_emu <- sqrt((10^etaEmu)*(k_ci^1.5))
    MSEs_Delta[holdout] <- mean((Delta_emu - Delta_holdout)^2)
    
    # plot(log10(k_ci), Delta_emu, type="l", 
    #      ylab = expression(Delta), main = paste("Delta", run))
    # lines(log10(k_ci), Delta_holdout, col="red")
    # legend("bottomright", legend = c("holdout","emu"), 
    #        lty=c(1,1), col = c("red","black"))
    if(holdout==1){
      plot(log10(k_ci), Delta_emu/Delta_holdout, type="l", lty=1, 
              xlab=expression(log[10](k)), ylab="Delta_emu/Delta_holdout",
              main = paste0("Interpolated CAMB runs, train on des=2"), 
              col=holdout, ylim=c(0.945,1.055))
      abline(h=1, col="black", lwd=2, lty=2)
    } else{
      lines(log10(k_ci), Delta_emu/Delta_holdout,
           col=holdout, ylim=c(0.945,1.055))
      print(Delta_emu/Delta_holdout)
    }
  } else {
    matplot(log10(k_ci), pk_cis, type="l", lty=1, 
            xlab=expression(log[10](k)), ylab="script P(k)",
            main = paste0("Interpolated CAMB runs, heldout=",holdout), col="gray")
    
    lines(log10(k_ci), pk_cis[,holdout], col="red")
    lines(log10(k_ci), etaEmu, col="blue", lty=2)
    legend("bottom", legend = c("runs","holdout","emu"), 
           lty=c(1,1,2), col = c("gray","red","blue"))
  }
}

# Omega_b is 4th column
des_test <- read.table("../CosmicEmu_etc/32 Model Runs CAMB/camb_new.design")
omega_b <- des_test[,4]
plot(omega_b, MSEs, main = "omega_b vs MSE (script P)")
plot(omega_b, MSEs_Delta, main = "omega_b vs MSE (Delta)")
if(PDF) dev.off()
