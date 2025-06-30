#############################################################
# compare predictions: DGP and GP-PC approach vs cosmic emu #
# also compare estimation method from process convolution.  #
#############################################################

library(dgp.hm)

# make JPEGs of plots?
JPG <- F

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
        col=cbcols[7], lwd=2, lty=3)
  # plot the cosmicEmu prediction
  lines(log10(k),cosmscrPsz - meanterm, col=cbcols[3], lwd=2, lty=2)
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

if(JPG) jpeg("../paper/pred_1to6.jpeg", width = 12, height = 4, units = "in", res = 300)
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
  meanterm = post_means_test[,model_name]#loess(y_avg ~ log10(k), span = 0.1)$fitted
  
  # calculate MSEs for each model
  mses[i,] <- (etaEmu[,model_name] - meanterm)^2
  mses_emu[i,] <- (cosmscrPsz - meanterm)^2
  mses_conv_t[i,] <- (as.numeric(convln_pred[i,]) - meanterm)^2
  mses_dgp_t[i,] <- (post_means_test[,model_name] - meanterm)^2

  # Make a plot showing each computer model, along with predictions
  # Average is removed to better visualize the variability and bias
  # par(mar=c(4,4.5,1,1), mfrow=c(1,1))
  # plot our dgp.fco prediction
  plot(log10(k), post_means_test[,model_name] - meanterm, col=cbcols[1], lwd=1,
       ylim = c(-.07,.07), type="l",
          xlab=expression(log[10](k)),
       yaxt=ifelse(i %in% c(2,3,5,6),"n","s"),
       xaxt=ifelse(i %in% c(1,2,3),"n","s"),
          ylab='script P')
  # plot the cosmic emu prediction
  lines(log10(k), cosmscrPsz - meanterm, col=cbcols[3], lty=3)
  # # plot the emu train
  # lines(log10(k), as.numeric(convln_pred[i,]) - meanterm, col=cbcols[6], lty=2, lwd=2)
  # plot the dgp.fco train
  lines(log10(k), etaEmu[,model_name] - meanterm, type="l", col = cbcols[7], lty=2)
  
  mtext(model_name, line = -2, adj=0.1)
  
}
mtext(expression(log[10](k)),side=1,line=2.75,outer=T)
mtext(paste0("Deviation from Posterior Mean"),side=2,line=2.5,outer=T)
if(JPG) dev.off()

# compare MSEs for each k value, by method
if(JPG) jpeg("../paper/mse_by_k.jpeg", width = 12, height = 4, units = "in", res = 300)
par(mfrow=c(1,1),oma=c(4,4,1.5,1),mar=c(0,0,0,0))
# plot cosmic emu col=3 (blue)
plot(log10(k),colMeans((mses_emu)), type="l", col=cbcols[3], lwd=2, lty=3,
     ylim = range(0,colMeans((mses)),colMeans((mses_emu))), ylab="")
# plot our dgp.fco col=7 (orange)
lines(log10(k),colMeans((mses)), col=cbcols[7], lwd=2, lty=2)
legend("topright",legend = c("Cosmic Emu","DGP.FCO"),
       col = c(cbcols[3],cbcols[7]), lty=c(3,2), lwd=c(2,2))
mtext(expression(log[10](k)),side=1,line=2.5,outer=T)
mtext(paste0("Average MSE"),side=2,line=2.5,outer=T)
if(JPG) dev.off()

# Make plot of basis function decomposition
# run predict.R code below first to load des_test and aps[[1]] 
# mean_PCs_oneW; hydro_plotting.R has original with GSMF
load("../pca/trained_GP_for_pca.rda")

# Read in test design
theta_fixed <- rep(0.5, 7)  # or use colMeans(des_test)[2:8] or something else
x1_seq <- seq(0, 1, length.out = 100)  # sweeping over theta_1
des_test <- cbind(x1_seq, matrix(rep(theta_fixed, each = 100), ncol = 7))
ntest <- nrow(des_test)

# make the predictions for the given des_test design
aps <- list()
for (i in 1:n_pc) {
  print(paste("GP pred",i))
  aps[[i]] = GPfit::predict.GP(as[[i]], des_test)
}

# create matrix of overall mean to add back on to predictions
mean_pred = matrix(mean_mat[,1],nrow=n_k,ncol=nrow(des_test))

# scale each PC by its predicted weight
eta_preds <- list()
for (i in 1:n_pc) eta_preds[[i]] = outer(bases[,i],aps[[i]]$Y_hat)

# combine all weighted PCs and the overall mean
etaEmu <- mean_pred + Reduce('+', eta_preds)

# mean_PCs_oneW.jpeg; PCA picture
if(JPG) jpeg("../paper/mean_PCs_oneW.jpeg", width = 12, height = 4, units = "in", res = 300)
par(mfrow=c(1,3), mar=c(4,4,1.8,1))
plot(log10(k), apply(eta,1,mean),type="l", xlab=expression(log[10](k)), ylab='script P')
matplot(log10(k), bases, type="l", xlab=expression(log[10](k)), ylab = 'deviation from average')
x1_pred <- des_test[,1]
plot(x1_pred[order(x1_pred)], aps[[5]]$Y_hat[order(x1_pred)], type="l", 
     xlab = expression(psi[1]), ylab = expression(gamma[5]))
if(JPG) dev.off()
