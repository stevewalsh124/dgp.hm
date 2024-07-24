# Power Spectra Simulations and linear theory (LT)

library(dgp.hm)

# pdf("LT_sims_scrP_all.pdf")

# All the Runs and Models are in units Mpc/h and h/Mpc.
# To be able to plot everything in one plot, we have to convert to Mpc or 1/Mpc.
# To convert k (column 1), we need to multiply it by hubble (to get k in units of 1/Mpc) 
# while the power spectrum (column 2) needs to be divided by h^3 (to get it in Mpc^3).

# contains model/run index, D+ values, and hubble values
data <- read.csv("../CosmicEmu_etc/32 Models Mira Titan/models_d+_h.csv")
Dp <- data$D.
h <- data$Hubble

# Working with runs 1-32 (aka models 11-42), but skip 17-18 (27-28)
run <- 1
model <- run + 10

# read in the linear output
sim <- read.table(paste0("../CosmicEmu_etc/32 Models Mira Titan/linear/RUN",
                         run,"/oneh_matterpower.dat"))
k_s <- sim$V1*h[run]
pk_s <- sim$V2/ (h[run])^3
plot(log10(k_s), scrP(pk_s, k_s), type="l", col="blue", lwd=4,
     xlab=expression(log[10](k)), ylab="script P(k)", ylim=c(-0.2,1.4), xlim = c(-2.3,0),
     main = paste0("LT and simulated spectra, model ", model))

# read in first batch of model runs (model #11)
pk_ms <- matrix(NA, 886, 15)
for (i in 0:15) {
  mod <- read.table(paste0("../CosmicEmu_etc/32 Models Mira Titan/pow.ic/M0", 
                           model, "/L1300/PM0",if(i < 10){"0"}, i,
                           "/analysis/Pow/m0", model, ".pk.ini"))
  k_m <- mod$V1 * h[model-10]
  pk_ms[,i] <- pk_m <- mod$V2 / (h[model-10])^3 / (Dp[model-10])^2
  lines(log10(k_m)[k_m <= 1], scrP(pk_m, k_m)[k_m <= 1], lty=3, col="gray")
}

lines(log10(k_m)[k_m <= 1], scrP(rowMeans(pk_ms), k_m)[k_m <= 1], 
      col="orange", lwd=2, lty=1)

legend(x = "bottom", legend = c("linear theory", "sims", "sim avg"),
       lty = c(1, 3, 1), lwd=c(4, 1, 2), col=c("blue","gray","orange"))


# Estimate the variance in a log-log regression model
# Should we estimate log-log regression using only the blue, or all red?

# only look at kvals between 0.01 and 1 (so, log10(k) is in [-2,0])
sub <- which(log10(k_m) <= 0 & log10(k_m) >= -2)
k_sub <- k_m[sub]

# obtain estimate of the model sims' variances at each k_sub value
var_pkms <- apply(pk_ms[sub,], 1, var)

# plot these log10(variances) against log10(k)
plot(log10(k_sub), log10(var_pkms),
     main = "log10(k) vs log10(Var(P(k)))")

# Get simple linear regression fit
lm_var <- lm(log10(var_pkms) ~ log10(k_sub))
lines(log10(k_sub), lm_var$fitted.values, col="lightblue", lwd=3)

# Get quadratic fit for log-log regression model
lm_var2 <- lm(log10(var_pkms) ~ log10(k_sub) + I(log10(k_sub)^2))
lines(log10(k_sub), lm_var2$fitted.values, col="red", lwd=3)

# Look at these variances on the original scale (not log10)
plot(k_sub, var_pkms, main = " k vs  Var(P(k))")
lines(k_sub, 10^(lm_var$fitted.values), col="lightblue", lwd=3)
lines(k_sub, 10^(lm_var2$fitted.values), col="red", lwd=3, lty=2)

# Add a legend
legend(x = "topright", 
       legend = c("estimated variances", "quad estimate", "SLR estimate"),
       lty = c(NA,2,1), pch=c(1,NA,NA), lwd=c(NA,2,2), col=c(1,2,"lightblue"))

# # Interpolate linear function onto the models' k values
# pk_si <- approx(k_s, scrP(pk_s, k_s), k_m)
# lines(log10(k_m), pk_si$y, col="red")

 # for (run in 2:32) {
#   model <- run + 10
#   if(model == 27 | model == 28) next
#   if(model < 20)  folder <- ""
#   if(model >= 20) folder <- "_2"
#   if(model >= 30) folder <- "_3"
#   if(model >= 40) folder <- "_4"
#   
#   sim <- read.table(paste0("../CosmicEmu_etc/32 Models Mira Titan/linear/RUN",
#                            run,"/oneh_matterpower.dat"))
#   k_s <- sim$V1*h[run]
#   pk_s <- sim$V2/ (h[run])^3
#   plot(log10(k_s), scrP(pk_s,k_s), type="l", col=run,
#        xlab=expression(log[10](k)), ylab="script P(k)", ylim=c(-0.2,1.4), xlim = c(-2.3,0),
#        main = paste0("LT and simulated spectra, model ", model))
#   
# 
#   for (i in 0:9) {
#     mod <- read.table(paste0("../CosmicEmu_etc/32 Models Mira Titan/pow",folder,
#                              ".ic/M0", model, 
#                              "/L1300/PM00",i,"/analysis/Pow/m0", model, ".pk.ini"))
#     k_m <- mod$V1 * h[model-10]
#     pk_m <- mod$V2 / (h[model-10])^3 / (Dp[model-10])^2
#     lines(log10(k_m)[k_m <= 1], scrP(pk_m, k_m)[k_m <= 1], lty=3, col=run)
#   }
# 
#   lines(log10(k_s), scrP(pk_s,k_s), col=run, lwd=2)
# 
# }

# dev.off()
