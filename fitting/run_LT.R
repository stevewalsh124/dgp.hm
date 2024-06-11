# Power Spectra Simulations and linear theory (LT)

library(dgp.hm)

# All the Runs and Models are in units Mpc/h and h/Mpc.
# To be able to plot everything in one plot, we have to convert to Mpc or 1/Mpc.
# To convert k (column 1), we need to multiply it by hubble (to get k in units of 1/Mpc) 
# while the power spectrum (column 2) needs to be divided by h^3 (to get it in Mpc^3).

# contains model/run index, D+ values, and hubble values
data <- read.csv("../CosmicEmu_etc/32 Models Mira Titan/models_d+_h.csv")
Dp <- data$D.
h <- data$Hubble

# From the "linear" folder
run <- 1
model <- run + 10
# read in first linear
sim <- read.table(paste0("../CosmicEmu_etc/32 Models Mira Titan/linear/RUN",
                         run,"/oneh_matterpower.dat"))
k_s <- sim$V1*h[run]
pk_s <- sim$V2/ (h[run])^3
# read in first model run
mod <- read.table(paste0("../CosmicEmu_etc/32 Models Mira Titan/pow.ic/M0", model, 
                         "/L1300/PM000/analysis/Pow/m0", model, ".pk.ini"))
k_m <- mod$V1 * h[model-10]
pk_m <- mod$V2 / (h[model-10])^3 / (Dp[model-10])^2
plot(log10(k_s), scrP(pk_s, k_s), type="l",
     xlab=expression(log[10](k)), ylab="script P(k)", ylim=c(-0.2,1.4), xlim = c(-2.3,0),
     main = "LT and simulated spectra")
lines(log10(k_m)[k_m <= 1], scrP(pk_m, k_m)[k_m <= 1], lty=3)

for (run in 2:32) {
  model <- run + 10
  if(model == 27 | model == 28) next
  if(model < 20)  folder <- ""
  if(model >= 20) folder <- "_2"
  if(model >= 30) folder <- "_3"
  if(model >= 40) folder <- "_4"
  mod <- read.table(paste0("../CosmicEmu_etc/32 Models Mira Titan/pow",folder,
                           ".ic/M0", model, 
                           "/L1300/PM000/analysis/Pow/m0", model, ".pk.ini"))
  k_m <- mod$V1 * h[model-10]
  pk_m <- mod$V2 / (h[model-10])^3 / (Dp[model-10])^2
  lines(log10(k_m)[k_m<=1], 
        scrP(pk_m,k_m)[k_m<=1], col=run, lwd=1, lty=3)
  
  sim <- read.table(paste0("../CosmicEmu_etc/32 Models Mira Titan/linear/RUN",
                           run,"/oneh_matterpower.dat"))
  k_s <- sim$V1*h[run]
  pk_s <- sim$V2/ (h[run])^3
  lines(log10(k_s), 
        scrP(pk_s,k_s), col=run, lwd=1)
}

legend(x = "bottom", legend = c("linear theory", "sims"),
       lty = c(1,3), lwd=c(1,1))


