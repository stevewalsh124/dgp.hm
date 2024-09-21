#############################################################
# compare predictions: DGP and GP-PC approach vs cosmic emu #
# also compare estimation method from process convolution.  #
#############################################################

library(dgp.hm)

# make a PDF of plots?
PDF <- F

# load the precision data (k, prec_highres, prec_lowres, index_list)
load("../Mira-Titan-IV-Data/precision_and_indexes.Rdata")

# pick which prediction scheme to look at
# number of principal components, and power for powered exponential cov fn
n_pc <- 10
pwdExp <- 1.9

args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))

# color-blind friendly color palette
cbcols <- palette.colors(palette = "Okabe-Ito")
#1 black, 2 orange, 3 skyblue , 4 bluishgreen, 5 yellow 
#6 blue, 7 vermillion,  8 reddishpurple, 9 gray 

if(PDF) pdf(paste0("compare_pred_",n_pc,"_",pwdExp,".pdf"))

# read in predicted power matter spectra for 6 test models
etaEmu <- read.csv(paste0("../pca/etaEmu",n_pc,"_",pwdExp,".csv"))
colnames(etaEmu) <- c("E001", "E002", "E003", "E009", "E010", "M000")

# import Kelly and Jared's process convolution estimates
# this is not emulating, this is estimating (so it should do better)
convln_pred <- read.csv("../CosmicEmu_etc/z_zero_test.txt", sep = "", header=F)

# record MSEs for each of four approaches
mses <- mses_emu <- mses_conv_t <- mses_dgp_t <- matrix(NA, 6, 351)

#####################################
# Plot predictions for both methods #
#####################################

par(mfrow=c(2,3),oma=c(4,4,1.5,1),mar=c(0,0,0,0))
# For each of the six test models...
for(i in 1:6){
  
  # Specify testing model name
  model <- c(112:116,0)[i]
  if (model <= 111) {
    model_name <- paste0("M", if (model < 100) {"0"}, if (model < 10) {"0"}, model) 
  } else {
    test_names <- c("E001", "E002", "E003", "E009", "E010")
    model_name <- test_names[model - 111]
  }
  
  # load test file's computer model runs
  # 1st column is k, 2nd is linear pert theory, 3:18 is low-res, 19 is hi-res
  file_name <- paste0("../Mira-Titan-IV-Data/Mira-Titan-2021/STEP499/pk_", 
                      model_name, "_test.dat")
  pk2 <- read.table(file_name)
  
  # load the weighted average from the computer model runs
  y_avg <- read.csv(paste0("../CosmicEmu_etc/y_avg_",model,".csv"))$x
  
  # load the test file's posterior mean (excluded from the training set)
  post_means_test <- read.csv("../fitting/results/post_means_test.csv")
  
  # read in Cosmic Emu prediction, then convert to script P
  cosmicEmu <- read.csv(paste0("../CosmicEmu_etc/EMU",
                               i-1,".txt"), sep="", header = F)
  cosmscrPsz <- scrP(cosmicEmu[,2],k)
  
  # select a mean term to make plots and calculate MSEs
  meanterm = y_avg
  
  # Make a plot showing each computer model, along with predictions
  # Average is removed to better visualize the variability and bias
  # par(mar=c(4,4.5,1,1), mfrow=c(1,1))
  matplot(log10(k), scrP(pk2[,3],k) - meanterm, type="n",
       ylim = c(-.07,.07),
       xlab=expression(log[10](k)),
       ylab='script P', yaxt=ifelse(i %in% c(2,3,5,6),"n","s"),
       xaxt=ifelse(i %in% c(1,2,3),"n","s"))
  # plot the low-res runs
  for (j in 3:18) lines(log10(k)[index_list$lowres.ix], 
                        (scrP(pk2[,j],k) - meanterm)[index_list$lowres.ix], 
                        col=cbcols[9], lwd=1)
  # plot the high res run
  lines(log10(k)[index_list$highres.ix], 
        (scrP(pk2[,19],k) - meanterm)[index_list$highres.ix], 
        col=cbcols[9], lwd=2)
  # plot the pert theory
  lines(log10(k)[index_list$pert.ix], 
        (scrP(pk2[,2],k) - meanterm)[index_list$pert.ix], 
        col=cbcols[9], lwd=3)
  # plot our prediction using DGPs and then GPs on weights of PCs (DGP-GPPC)
  lines(log10(k), (etaEmu[,model_name]) - meanterm, 
        col=cbcols[3], lwd=2, lty=3)
  # plot the cosmicEmu prediction
  lines(log10(k),cosmscrPsz - meanterm, col=cbcols[7], lwd=2, lty=2)
  # # plot the proc convln estimate (trained with model runs)
  # lines(log10(k), convln_pred[which(c(112:116,0)==model),] - meanterm, 
  #       col=cbcols[3],lty=2, lwd=2)
  # # plot the dgp.hm estimate (trained with model runs)
  # lines(log10(k), post_means_test[,model_name] - meanterm, 
  #       col=cbcols[6],lty=2, lwd=2)
  # Make a legend
  # legend("topright",legend = c("data","cosmicEMU","DGP-GPPC"),
  #        col = c(cbcols[9],cbcols[7],cbcols[3]), 
  #        lty=c(1,1,1,2,3,2), lwd=2)
  mtext(model_name, line = -2, adj=0.1)
  
  # calculate MSEs for each model
  mses[i,] <- (etaEmu[,model_name] - meanterm)^2
  mses_emu[i,] <- (cosmscrPsz - meanterm)^2
  mses_conv_t[i,] <- (as.numeric(convln_pred[i,]) - meanterm)^2
  mses_dgp_t[i,] <- (post_means_test[,model_name] - meanterm)^2
  # cat(c("RMSE_dgp",sqrt(sum((etaEmu[,i] - meanterm)^2))),"\n")
  # cat(c("RMSE_emu",(sqrt(sum((cosmscrPsz - meanterm)^2)))),"\n")
  # cat(c("RMSE_conv_t",sqrt(sum((convln_pred[i,] - meanterm)^2))),"\n")
  # cat(c("RMSE_dgp_t",sqrt(sum((post_means_test[,model_name] - meanterm)^2))),"\n")
  
}
mtext(expression(log[10](k)),side=1,line=2.75,outer=T)
mtext(paste0("Average RMSE"),side=2,line=2.5,outer=T)

# compare RMSEs for each k value, by method
par(mfrow=c(1,2),oma=c(4,4,1.5,1),mar=c(0,0,0,0))

plot(log10(k),colMeans(sqrt(mses_emu)), type="l", col=cbcols[7], lwd=2,
     ylim = range(0,colMeans(sqrt(mses)),colMeans(sqrt(mses_emu))))
lines(log10(k),colMeans(sqrt(mses)), col=cbcols[3], lwd=2)
legend("topright",legend = c("cosmicEMU","DGP.HM"),
       col = c(cbcols[7],cbcols[3]), lwd=c(2,2))
plot(log10(k),colMeans(sqrt(mses_conv_t)), col=cbcols[7], type="l",
     yaxt="n", ylim = range(0,colMeans(sqrt(mses))))
lines(log10(k),colMeans(sqrt(mses_dgp_t)), col=cbcols[3])
legend("topright",legend = c("DPC","DGP.HM"),
       col = c(cbcols[7],cbcols[3]), lwd=c(1,1))

mtext(expression(log[10](k)),side=1,line=2.5,outer=T)
mtext(paste0("Average RMSE"),side=2,line=2.5,outer=T)


# boxplot of RMSEs by method
par(mfrow=c(1,1))
boxplot(c(sqrt(mses)),c(sqrt(mses_emu)),c(sqrt(mses_conv_t)),c(sqrt(mses_dgp_t)), 
        names=c("DGP-GPPC","EMU","conv_t","dgp_t"), xlab="RMSE", horizontal = T,
        col = c(cbcols[5],cbcols[2],cbcols[3],cbcols[6]))

# compare average MSE for both prediction methods
mean(mses_emu)
mean(mses)

# find proportion of k values where our method has lower MSE
mean(colMeans(mses) < colMeans(mses_emu))

# Make plot of basis function decomposition
# run predict.R first to load des_test and aps[[1]]
# mean_PCs_oneW; hydro_plotting.R has original with GSMF
par(mfrow=c(1,3), mar=c(4,4,1.8,1))
plot(apply(eta,1,mean),type="l", xlab=expression(log[10](k)), ylab='script P')
matplot(bases, type="l", xlab=expression(log[10](k)), ylab = 'deviation from average')
des_test <- read.csv("../Mira-Titan-IV-Data/design_test.txt",sep="",header = F)
x1_pred <- des_test[,1]
plot(x1_pred[order(x1_pred)], aps[[1]]$Y_hat[order(x1_pred)], type="l", 
     xlab = expression(theta[1]), ylab = expression(W[1]))

if(PDF) dev.off()
