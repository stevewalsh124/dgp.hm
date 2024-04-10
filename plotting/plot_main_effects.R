#####################################
# Script with plotting for PCA fits #
# and main effects estimates, and   #
# variance decomposition            #
#####################################

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

# Read in test design
# des_test <- read.csv("../Mira-Titan-IV-Data/design_test.txt",sep="",header = F)
m = 100
des_test <- read.csv(paste0("../csv/des_test_unif",m,".csv"), row.names = 1)
ntest <- nrow(des_test)

# read in predicted power matter spectra for test models
etaEmu <- read.csv(paste0("../pca/etaEmu_unif",m,".csv"))
# colnames(etaEmu) <- c("E001", "E002", "E003", "E009", "E010", "M000")
ntest = ncol(etaEmu)

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

# number of low-to-high variations for each input
# to visualize the main effects
nv = 9 
prange = seq(0,1,length=nv)

# plot them
matplot(kvals, eta, type="l", 
        xlab = "log10(Stellar Mass)", ylab = "log10(GSMF_Apperture)", 
        main = paste0(nruns, " training runs"))

plot(a$d)   # looks like 4+ pc's should work

# plot the basis elements
matplot(kvals,a$u%*%sqrt(diag(a$d)),type='l')

hist(coef1)

matplot(bases, type="l")

par(mfrow=c(2,2),oma=c(0,0,0,0),mar=c(4,4,1.5,1))
# use all 111 bases
matplot(kvals,bases%*%t(coef),type='l',ylab='P(k)')
mtext(paste0('all (',ncol(bases),') bases'),side=3,lin=.2)
# use n_pc bases
matplot(kvals,bases[,1:n_pc]%*%t(coef[,1:n_pc]),type='l',ylab='P(k)')
mtext(paste0(n_pc,' bases'),side=3,lin=.2)

# note, the spectra can be recovered by adding back in the mean
matplot(kvals,bases%*%t(coef)+mean_mat,type='l',ylab='P(k)')
mtext(paste0('all (',ncol(bases),') bases'),side=3,lin=.2)
# use n_pc bases
matplot(kvals,bases[,1:n_pc]%*%t(coef[,1:n_pc])+mean_mat,type='l',ylab='P(k)')
mtext(paste0(n_pc,' bases'),side=3,lin=.2)

# show emulator output compared to training data
par(mfrow=c(1,2),oma=c(0,0,0,0),mar=c(4,4,1.5,1))
matplot(kvals,eta,type='l',ylab='GSMF_A')
mtext(paste0(nruns,'-run training set'),side=3,line=.1,cex=.9)

matplot(kvals,etaEmu,type='l',ylab='emulator')
mtext(paste0(ntest,'-run design, GP-PC'),side=3,line=.1,cex=.9)

# compute the main effects for each param
mainEffs_gppc = array(NA,c(n_k,nv,p))
for(k in 1:p){
  # collect GP-PC predictions
  sims1 = as.matrix(etaEmu[,(k-1)*(ntest/p)+1:(ntest/p)])
  # reshape into a 3d array by kval,rep, paramval
  sims1a = array(sims1,c(n_k,m,nv))
  # compute main effect
  mainEffs_gppc[,,k] = apply(sims1a,c(1,3),mean)
}

avg_mean <- rowMeans(eta)

# make a plot for the main effects
# set the parameter names
# c('omega_m','omega_b','sigma_8','h','n_s','w_0','w_a','omega_nu','z')
pnames = c(paste0("param",1:p))

par(mfrow=c(2,4),oma=c(4,4,1.5,1),mar=c(0,0,0,0))

# make a grayscale
grcolors = paste0("grey",round(seq(90,25,length=11)),sep='')

# Plot GP-PC main effects (average removed)
yr = range(mainEffs_gppc - avg_mean)
for(k in 1:p){
  matplot(kvals,mainEffs_gppc[,,k] - avg_mean,
          type='l',col=grcolors,lty=1,ylim=yr,axes=F); box()
  text(-3,yr[1]+0.01,pnames[k],adj=c(0,0))
  if(k %in% c(1,5)) axis(2)
  if(k %in% c(5:p)) axis(1)
}
mtext('log10(Stellar Mass)',side=1,line=2.5,outer=T)
mtext(paste0("log10(GSMF_A) [average removed top row]"),side=2,line=2.5,outer=T)
mtext("GP-PC main effect estimates", side=3,outer = T)

# Plot GP_PC main effects (average not removed)
yr = range(mainEffs_gppc)
for(k in 1:p){
  matplot(kvals,mainEffs_gppc[,,k],
          type='l',col=grcolors,lty=1,ylim=yr,axes=F); box()
  text(-3,yr[1]+0.01,pnames[k],adj=c(0,0))
  if(k %in% c(1,5)) axis(2)
  if(k %in% c(5:p)) axis(1)
}
mtext('log10(Stellar Mass)',side=1,line=2.5,outer=T)
mtext(paste0('log10(GSMF_A)'),side=2,line=2.5,outer=T)
mtext("GP-PC main effect estimates", side=3,outer = T)

# expand the test design to incorporate k
des_testK <- matrix(NA, nrow = ntest*n_k, ncol = p+1)
for (i in 1:ntest) {
  for (j in 1:n_k) {
    des_testK[(i-1)*(n_k) + j,] <- as.numeric(cbind(kvals[j], des_test[i,]))
  }
}

# appropriately reshape the model output
yEmu = as.vector(as.matrix(etaEmu))

# convert columns of des_testK to factors
for(k in 1:ncol(des_testK)) des_testK[,k] = as.factor(des_testK[,k])
colnames(des_testK) = c("k",paste0("p",1:p))

# run the linear model to estimate effects
# everything interacts with k.
facSAdf = data.frame(y1=yEmu, des_testK)
aEmu = lm(y1 ~ k+k:(p1+p2+p3+p4+p5+p6+p7+p8)^3,data=facSAdf)

# get the Sum of Squares for each term using the anova() function
bEmu = anova(aEmu)

# make a plot of the sums of squares
ssnames = row.names(bEmu)
ss = ssEmu = bEmu$`Sum Sq`

ss_lim <- range(sqrt(c(ss)))

# plot sum of squares (ss) decomposition
par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(4,4,1.8,1))
plot(1:length(ss), sqrt(ssEmu), xlab='effect', ylab='emu sqrt SS',
     col='blue', log='y', ylim=ss_lim)
mtext('GP-PC emulator', side=3, line=.2, cex=.9)
abline(v=c(1.5, 9.5, 37.5, 93.5), col="gray", lty=2)
abline(h=1, col="red", lty=2)

# identify large effects on the plot
# plot(1:2^p,sqrt(ssEmu),xlab='effect',ylab='emu+error sqrt SS',col='red',log='y')
# identify(1:2^p,sqrt(ssEmu), labels = ssnames)
# ssnames[c(1,2,11,16)] # residuals is important
