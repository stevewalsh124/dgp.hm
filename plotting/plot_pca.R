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

if(PDF) pdf(paste0("compare_pred_",n_pc,"_",pwdExp,".pdf"))

# read in predicted power matter spectra for 6 test models
etaEmu <- read.csv(paste0("../pca/etaEmu",n_pc,"_",pwdExp,".csv"))
colnames(etaEmu) <- c("E001", "E002", "E003", "E009", "E010", "M000")

# import Kelly and Jared's process convolution estimates
# this is not emulating, this is estimating (so it should do better)
convln_pred <- read.csv("../CosmicEmu_etc/z_zero_test.txt", sep = "", header=F)

# record MSEs for each of three approaches
mses <- mses_emu <- mses_conv_t <- mses_dgp_t <- matrix(NA, 6, 351)

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
  par(mar=c(4,4.5,1,1), mfrow=c(1,1))
  matplot(log10(k), scrP(pk2[,3],k) - meanterm, type="n",
       ylim = c(-.12,.12),
       xlab=expression(log[10](k)),
       ylab='script P')
  # plot the low-res runs
  for (j in 3:18) lines(log10(k)[index_list$lowres.ix], 
                        (scrP(pk2[,j],k) - meanterm)[index_list$lowres.ix], 
                        col="gray", lwd=1)
  # plot the high res run
  lines(log10(k)[index_list$highres.ix], 
        (scrP(pk2[,19],k) - meanterm)[index_list$highres.ix], 
        col="red", lwd=2)
  # plot the pert theory
  lines(log10(k)[index_list$pert.ix], 
        (scrP(pk2[,2],k) - meanterm)[index_list$pert.ix], 
        col="green", lwd=3)
  # plot our prediction using DGPs and then GPs on weights of PCs (DGP-GPPC)
  lines(log10(k), (etaEmu[,model_name]) - meanterm, 
        col="magenta", lwd=2, lty=3)
  # plot the cosmicEmu prediction
  lines(log10(k),cosmscrPsz - meanterm, col="orange", lwd=2, lty=2)
  # plot the proc convln estimate (trained with model runs)
  lines(log10(k), convln_pred[which(c(112:116,0)==model),] - meanterm, 
        col="pink",lty=2, lwd=2)
  # plot the proc convln estimate (trained with model runs)
  lines(log10(k), post_means_test[,model_name] - meanterm, 
        col="blue",lty=2, lwd=2)
  # Make a legend
  legend("topright",legend = c("low res","pert","hi res","cosmicEMU","DGP-GPPC","conv_t","dgp_t"),
         col = c("gray","green","red","orange","magenta","pink","blue"), 
         lty=c(1,1,1,2,3,2), lwd=2)
  mtext(model_name)

  # calculate MSEs for each model
  mses[i,] <- (etaEmu[,model_name] - meanterm)^2
  mses_emu[i,] <- (cosmscrPsz - meanterm)^2
  mses_conv_t[i,] <- (as.numeric(convln_pred[i,]) - meanterm)^2
  mses_dgp_t[i,] <- (post_means_test[,model_name] - meanterm)^2
  cat(c("RMSE_dgp",sqrt(sum((etaEmu[,i] - meanterm)^2))),"\n")
  cat(c("RMSE_emu",(sqrt(sum((cosmscrPsz - meanterm)^2)))),"\n")
  cat(c("RMSE_conv_t",sqrt(sum((convln_pred[i,] - meanterm)^2))),"\n")
  cat(c("RMSE_dgp_t",sqrt(sum((post_means_test[,model_name] - meanterm)^2))),"\n")
  
}

# compare RMSEs for each k value, by method
plot(log10(k),colMeans(sqrt(mses)), type="l", col="magenta", lwd=2,
     ylim = range(0,colMeans(sqrt(mses)),colMeans(sqrt(mses_emu))),
     main = "average RMSE by k value")
lines(log10(k),colMeans(sqrt(mses_emu)), col="orange", lwd=2)
lines(log10(k),colMeans(sqrt(mses_conv_t)), col="pink")
lines(log10(k),colMeans(sqrt(mses_dgp_t)), col="blue")
legend("topright",legend = c("cosmicEMU","DGP-GPPC","conv_t","dgp_t"),
       col = c("orange","magenta","pink","blue"), lwd=c(2,2,1,1))

# boxplot of RMSEs by method
boxplot(c(sqrt(mses)),c(sqrt(mses_emu)),c(sqrt(mses_conv_t)),c(sqrt(mses_dgp_t)), 
        names=c("DGP-GPPC","EMU","conv_t","dgp_t"), ylab="RMSE",
        main = "all RMSEs across k values and cosmologies", 
        col = c("magenta","orange","pink","blue"))

# compare average MSE for both prediction methods
mean(mses_emu)
mean(mses)

# find proportion of k values where our method has lower MSE
mean(colMeans(mses) < colMeans(mses_emu))

if(PDF) dev.off()
