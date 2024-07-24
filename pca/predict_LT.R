#######################################################################
# This script obtains principal components and their corresponding weights 
# Then predict cosmologies' spectra for LINEAR THEORY test cases
#
# Models 11-42 (except 27-28) are considered here.
#
###############################################################################

pdf("predict_LT_9PCs.pdf")

tic <- proc.time()[3]

library(GPfit)

args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))

#########################################
# Preprocess the data (h, interp, etc.) #
#########################################

# load the training design
all_runs <- c(1:16,19:32)
for(test_run in c(5,10,15,20,25,30)){ #1:16,19:32
  print(paste0("test run is ", test_run))
  training_runs <- setdiff(all_runs, test_run)
  des_train_full <- read.csv("../Mira-Titan-IV-Data/design_train.txt", header=F, sep = "")
  des_train <- des_train_full[training_runs,]
  p = ncol(des_train)  # number of input parameters
  
  # contains model/run index, D+ values, and hubble values
  # in order to appropriate scale the linear theory (and model runs)
  data <- read.csv("../CosmicEmu_etc/32 Models Mira Titan/models_d+_h.csv")
  Dp <- data$D.
  h <- data$Hubble
  
  # Load the posterior means obtained from run_fit.R and collect_means.R
  # Test on run #19
  k_int <- seq(-4, 1, by=0.01)
  eta <- matrix(NA, length(k_int), length(training_runs))
  for (i in 1:29) {
    run <- training_runs[i]
    sim <- read.table(paste0("../CosmicEmu_etc/32 Models Mira Titan/linear/RUN",
                             run,"/oneh_matterpower.dat"))
    k_s <- sim$V1*h[run]
    pk_s <- sim$V2/ (h[run])^3
    eta[,i] <- approx(log10(k_s), scrP(pk_s,k_s), k_int)$y 
  }
  matplot(k_int, eta, type="l", lty=1, col="gray", main = paste0("test #", test_run))
  nruns = ncol(eta)
  
  # Use log10 k values
  kvals <- k_int
  
  #################################
  # Obtain weights for PCs via GP #
  #################################
  
  # n_pc: number of principal components to model
  n_pc = 9
  if(n_pc < 1 | n_pc > 111) stop("number of PCs has to be between 1 and 111")
  
  # choose the power for powered-exponential covariance function 
  pwdExp = 1.9
  
  # n_k: number of k values considered for each cosmology
  n_k <- length(kvals)
  
  # calculate overall mean of 111 posterior means
  mean_mat = matrix(apply(eta,1,mean),nrow=n_k,ncol=nruns)
  
  # convert all 111 posterior means to have zero mean
  eta0 = eta - mean_mat
  
  # use singular value decomposition
  a = svd(eta0)
  
  # look at coefficients for each basis function
  coef1 = a$v[,1]
  
  # scale the coefficients so they have variance = 1
  coef = a$v*sqrt(nruns)
  
  # accordingly, scale the bases so when multiplied by
  # the coef's, we'll get the spectra back
  bases = a$u%*%diag(a$d)/sqrt(nruns)
  spectraFull = bases%*%t(coef)
  
  # try fitting gp's to the coefficients
  as <- list()
  for (i in 1:n_pc) {
    print(paste("GP",i))
    as[[i]] = GP_fit(des_train, coef[,i], 
                     corr = list(type="exponential", power=pwdExp))
  }
  
  #############################
  # Prediction of held-out LT #
  #############################
  
  # Read in test design
  des_test <- des_train_full[test_run,]
  ntest <- nrow(des_test)
  
  # make the predictions for the given des_test design
  aps <- list()
  for (i in 1:n_pc) {
    print(paste("GP pred",i))
    aps[[i]] = predict(as[[i]],des_test)
  }
  
  # create matrix of overall mean to add back on to predictions
  mean_pred = matrix(mean_mat[,1],nrow=n_k,ncol=nrow(des_test))
  
  # scale each PC by its predicted weight
  eta_preds <- list()
  for (i in 1:n_pc) eta_preds[[i]] = outer(bases[,i],aps[[i]]$Y_hat)
  
  # combine all weighted PCs and the overall mean
  etaEmu <- mean_pred + Reduce('+', eta_preds)
  lines(k_int, etaEmu)
  
  # write the predicted spectra
  # write.csv(etaEmu, paste0("etaEmu",n_pc,"_",pwdExp,"_LT.csv"), row.names = FALSE)
  
  # print time for training and predictions
  toc <- proc.time()[3]
  toc - tic
  
  ####################
  # Plot the results #
  ####################
  
  # plot the held out test run to compare with prediction
  for (run in test_run) {
    sim <- read.table(paste0("../CosmicEmu_etc/32 Models Mira Titan/linear/RUN",
                             run,"/oneh_matterpower.dat"))
    k_s <- sim$V1*h[run]
    pk_s <- sim$V2/ (h[run])^3
    lines(log10(k_s), scrP(pk_s, k_s), col="red", lty=2)
  }
  
  legend(x = "bottom", legend = c("training", "prediction", "held-out/truth"),
         col=c("gray","black","red"), lty = c(1,1,2), lwd=c(1,1,1))
  
  par(mfrow=c(1,3))
  betas <- c()
  for (i in 1:n_pc) betas <- c(betas, as[[i]]$beta)
  hist(betas, main = paste0("hist of betas for run ",test_run))
  hist(10^betas, main = paste0("10^betas for run ",test_run))
  hist(1/(10^betas), main = paste0("1/(10^betas) for run ",test_run))
  par(mfrow=c(1,1))
}

dev.off()