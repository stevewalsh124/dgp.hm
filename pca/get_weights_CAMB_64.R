###############################################################################
# This script visualizes and interpolates the CAMB model runs, and also 
# obtains principal components and their corresponding weights 
# These will then get used to predict spectra for holdout cases
#
# Models 1-32 are used here.
#
###############################################################################

tic <- proc.time()[3]

library(GPfit)
library(dgp.hm)

# pdf("CAMB_sims_scrP_all.pdf")

# All the Runs and Models are in units Mpc/h and h/Mpc.
# To be able to plot everything in one plot, we have to convert to Mpc or 1/Mpc.
# To convert k (column 1), we need to multiply it by hubble (to get k in units of 1/Mpc) 
# while the power spectrum (column 2) needs to be divided by h^3 (to get it in Mpc^3).

# load the training design (ignore 7th column, indicating des=1 or des=2)
X = read.csv("../CosmicEmu_etc/pow64/cambDesigns_32x6x2.csv")

# train on design 1, test on design 2
i_train = which(X$des == 2)
nruns = length(i_train)
des_train = X[i_train, -7]

p = ncol(des_train)  # number of input parameters
# hubble values to adjust k and P(k)
h <- read.csv("../CosmicEmu_etc/pow64/cambDesigns_32x6x2.csv")$h

# read in the linear output
run <- i_train[1]
sim <- read.table(paste0("../CosmicEmu_etc/pow64/RUN",
                         run,"/oneh_matterpower.dat"))
k_c <- sim$V1 * h[run]
pk_c <- sim$V2/ (h[run])^3
plot(log10(k_c), scrP(pk_c, k_c), type="l",
     xlab=expression(log[10](k)), ylab="script P(k)",
     xlim = c(-4.2 , 1.2), ylim=c(-5.2 ,1.6), 
     main = "held-out CAMB runs (adjusted by h)")

# create k values of CAMB for interpolation
k_ci <- 10^(seq(-1.7, 0.7, by=0.015))

# interpolate to common k grid (k_ci)
pk_ci <- approx(k_c, scrP(pk_c, k_c), k_ci)$y

# Collect training runs
# So, dim of pk_cis is (# of interp k vals) x (# of CAMB runs)
pk_cs <- matrix(NA, length(pk_c), nruns)
pk_cis <- matrix(NA, length(pk_ci), nruns)
k_cs <- matrix(NA, length(k_c), nruns)
pk_cis[,1] <- pk_ci
pk_cs[,1] <- pk_c
k_cs[,1] <- k_c

for (i in 2:nruns) {
  run = i_train[i]
  sim <- read.table(paste0("../CosmicEmu_etc/pow64/RUN",
                           run,"/oneh_matterpower.dat"))
  k_cs[,i] <- k_c <- sim$V1 * h[run]
  pk_cs[,i] <- pk_c <- sim$V2 / (h[run])^3
  lines(log10(k_c), scrP(pk_c, k_c), type="l", col=run)
  pk_cis[,i] <- approx(k_c, scrP(pk_c, k_c), k_ci)$y
  
}

# Plot the interpolated CAMB runs
matplot(log10(k_ci), pk_cis, type="l", lty=1, 
        xlab=expression(log[10](k)), ylab="script P(k)",
        main = "Interpolated CAMB runs")

# load the wavenumber vals
kvals <- log10(k_ci)

# Load the posterior means obtained from run_fit.R and collect_means.R
eta <- pk_cis
if(nruns != ncol(eta)) stop("nruns should equal ncol(eta)")

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

save.image(paste0("trained_GP_for_pca_CAMB",
                  "des1",".rda"))

