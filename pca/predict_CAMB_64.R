
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
i_test = which(X$des != 1)
Xp <- read.csv("../CosmicEmu_etc/pow64/cambDesigns_32x6x2.csv")[i_test,-7]

# standardize designs to be in the unit hypercube [0,1]^6
# ignore the last column of X, which denotes des=1 or des=2
for (i in 1:(ncol(X)-1)) Xp[,i] = (Xp[,i] - min(Xp[,i]))/(max(Xp[,i])-min(Xp[,i]))
ntest <- nrow(Xp)

# make the predictions for the given Xp design
aps <- list()
for (i in 1:n_pc) {
  print(paste("GP pred",i))
  aps[[i]] = predict(as[[i]],Xp)
}

# create matrix of overall mean to add back on to predictions
mean_pred = matrix(mean_mat[,1],nrow=n_k,ncol=nrow(Xp))

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
