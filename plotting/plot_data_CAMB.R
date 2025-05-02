
###############################################################################
# This script fits the Bayesian hierarchical model for a particular CAMB model
# and writes the results to a csv file.
#
# Command line arguments:
#    model: number of particular model (1-32)
#    deep: indicator for whether to fit a DGP or GP (1 = deep, 0 = not deep)
# 
# Outputs:
#    writes csv file of predicted values to the results folder
#
###############################################################################

library(dgp.hm)
library(zoo) # rollmean
library(Matrix) # bdiag
library(pracma) #interp1

load(file = "../fitting/results/CAMB/fit_1camb_image.rda")

# # Plot the first batch of model runs for CAMB
# plot(x_lo, y_lo[,1],
#      col="gray", type="l", main = model, xlim=c(-3,0), ylim=c(-2,2),
#      ylab = "script P (adj by h,g)",
#      xlab = "log10(k) (adj by h)")
# for (i in 2:n_lo) {
#   if(model==8 & i==9) next
#   lr = read.table(paste0("../CosmicEmu_etc/pows/M0",if(model<10){"0"},model,
#                          "/L1300/PM0",if(i<10){"0"},i,"/output/m0",
#                          if(model<10){"0"},model,".pk.ini"))
#   y_lo[,i] = scrP(lr[which(h*lr[,1]<=x_ub),2]/h^3/g^2, 10^x_lo)
#   lines(x_lo, y_lo[,i], col="gray")
# }
# lines(x_hi, y_hi, col="red")#, lwd=2)
# lines(x_camb, y_camb, col="blue1",type="l", lwd=3, lty=3)
# legend("bottomright", c("camb","lr","hr"), lty=1, col=c("blue","grey","red"))


# color-blind friendly color palette
cbcols <- palette.colors(palette = "Okabe-Ito")
#1 black, 2 orange, 3 skyblue , 4 bluishgreen, 5 yellow 
#6 blue, 7 vermillion,  8 reddishpurple, 9 gray 

# Plot showing dgp.hm fit compared with camb
# Plot with all of [-2.5, -2.2] camb
# CAMB_fit_model1
png("../paper/CAMB_fit_model1.png", width = 11, height = 4.5, 
    units = "in", res=300)
par(mar=c(4.5, 4.5, 1.5, 1.5) + 0.1)
plot(x, y_cambi - loess_fit$fitted, type="l", col=cbcols[3],
     ylim = range(y_cambi - loess_fit$fitted)+c(-.025,.025),
     xlab=expression(paste(log[10](k),", transformed to [0,1]")), 
     ylab=expression(paste("\U1D4AB","  (k), centered")), cex.lab = 1.35)
     # main = paste(model, "cover", coverages[model]))
abline(h=0, col="black", lty=2)
for (i in 1:n_lo) lines(x, y_loi[,i] - loess_fit$fitted, col=cbcols[9], lwd=0.4)
lines(x, y_hii - loess_fit$fitted, col=cbcols[8], lwd=1.5)
lines(x, y_avg - loess_fit$fitted, col=cbcols[4], lwd=2, lty=2)
lines(x, fit$m - loess_fit$fitted , col=cbcols[2], lwd=2)
lines(x, fit$ub - loess_fit$fitted, col=cbcols[2], lwd=2, lty=3)
lines(x, fit$lb - loess_fit$fitted, col=cbcols[2], lwd=2, lty=3)
lines(x, y_cambi - loess_fit$fitted, col=cbcols[3], lwd=2)
legend("topright", c("inf-res","low-res","hi-res","wt avg","dgp.hm"), 
       lty=c(1,1,1,2,1), cex=1, bty="n",
       col=c(cbcols[3],cbcols[9],cbcols[8],cbcols[4],cbcols[2]), lwd=c(2,0.5,1.5,2,2))
dev.off()
