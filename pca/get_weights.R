
###############################################################################
# This script obtains principal components and their corresponding weights 
# These will then get used to predict cosmologies' spectra for test cases
#
# Models 1-111 are saved in "post_means_train.csv".
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

#################################
# Obtain weights for PCs via GP #
#################################

# n_pc: number of principal components to model
n_pc = 10
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

# print time for training GPs for the weights
toc <- proc.time()[3]
toc - tic

save.image("trained_GP_for_pca.rda")