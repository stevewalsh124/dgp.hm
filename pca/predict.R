
###############################################################################
# This uses the PCs and weights (from get_weight.R) to predict cosmologies' 
# spectra for unobserved testing designs (6 hold out cosmologies)
#
# Models 1-111 are saved in "post_means_train.csv".
# Models 0,112-116 are saved in "post_means_test.csv".
#
###############################################################################

library(GPfit)

# Load in trained GPs for the weights of the PCs, and other necessary pieces
load("trained_GP_for_pca.rda")

tic <- proc.time()[3]

# Read in test design
des_test <- read.csv("../Mira-Titan-IV-Data/design_test.txt",sep="",header = F)
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

# write the predicted spectra
write.csv(etaEmu, paste0("etaEmu",n_pc,"_",pwdExp,".csv"), row.names = FALSE)

# print time for predictions
toc <- proc.time()[3]
toc - tic
