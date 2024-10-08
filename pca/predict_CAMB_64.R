
###############################################################################
# This uses the PCs & weights (from get_weights_CAMB.R) to predict cosmologies' 
# spectra for unobserved testing designs (1 hold out cosmologies)
#
# CAMB models 1-32 are used in a leave-one-out manner.
#
###############################################################################

library(GPfit)

# Load in trained GPs for the weights of the PCs, and other necessary pieces
load(paste0("trained_GP_for_pca_CAMB","des1",".rda"))

tic <- proc.time()[3]

# Read in test design
i_test = which(X$des != 2)
des_test <- read.csv("../CosmicEmu_etc/pow64/cambDesigns_32x6x2.csv")[i_test,-7]
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
write.csv(etaEmu, paste0("etaEmu",n_pc,"_",pwdExp,"CAMB","des1",".csv"), row.names = FALSE)

# print time for predictions
toc <- proc.time()[3]
toc - tic
