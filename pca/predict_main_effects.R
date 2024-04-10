
###############################################################################
# This script obtains batch predictions for each of the 6 test cosmologies
#
# Models 1-111 are saved in "post_means_train.csv".
# Models 0, 112-116 are saved in "post_means_test.csv".
#
###############################################################################

tic <- proc.time()[3]

library(GPfit)

# redshift of z
step <- 499 

args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))

##############################
# Load the preprocessed data #
##############################

# load the training design
des_train <- read.csv("../Mira-Titan-IV-Data/design_train.txt", header=F, sep = "")
p = ncol(des_train)  # number of input parameters

# load the wavenumber vals
load("../Mira-Titan-IV-Data/precision_and_indexes.Rdata")
kvals <- log10(k)

# Load the posterior means obtained from run_fit.R and collect_means.R
eta <- read.csv("../fitting/results/post_means_train.csv")
nruns = ncol(eta)

######################
# Obtain predictions #
######################

# n_pc: number of principal components to model
n_pc = 5
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
  as[[i]] = GP_fit(des_train,coef[,i], 
                   corr = list(type="exponential",power=pwdExp))
}

# Read in test design
m = 100
nv = 9
prange <- seq(0,1,len=nv)
# generate parameter settings that are uniform over [0,1]^8
udesign = matrix(runif(p*m),ncol=p)

# estimate main effect for parameter 1
# repeat the udesign matrix nv times
repUdes = matrix(rep(t(udesign),nv),ncol=p,byrow=TRUE)
# set up a vector that varies parameter value with each copy of udesign
repPrange = matrix(prange,ncol=m,nrow=nv)
repPrange = as.vector(t(repPrange))

# do this for all 8 parameters
list_facDes <- list()
for(k in 1:p){
  des1 = repUdes
  des1[,k] = repPrange
  colnames(des1) <- paste0("Var", 1:p)
  list_facDes[[k]] <- des1
}

facDes_all <- rbind(list_facDes[[1]],list_facDes[[2]],list_facDes[[3]],list_facDes[[4]],
                    list_facDes[[5]],list_facDes[[6]],list_facDes[[7]],list_facDes[[8]])
facDes_all <- unname(facDes_all)
des_test <- facDes_all

#read.csv("../Mira-Titan-IV-Data/design_test.txt",sep="",header = F)
write.csv(des_test, file=paste0("../csv/des_test_unif",m,".csv"))
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
write.csv(etaEmu, paste0("../pca/etaEmu_unif",m,".csv"), row.names = FALSE)

# print time for predictions
toc <- proc.time()[3]
toc - tic
