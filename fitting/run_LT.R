# Power Spectra Simulations and linear theory (LT)

library(dgp.hm)

# D+ values to scale sims in order to match the linear theory
Dp <- c(0.007985, 0.006853, 0.007923, 0.007024, 0.007214, 
        0.008146, 0.006776, 0.007998, 0.007517, 0.006175) 
        # 0.010293, 0.008178, 0.011243, 0.007489, 0.009507, 
        # 0.013198, 0.007817, 0.008683, 0.006918, 0.008124, 
        # 0.006601, 0.006948, 0.009904, 0.006510, 0.007727, 
        # 0.010022, 0.007219, 0.007437, 0.010621)

# From the "linear" folder
sim <- read.table("../CosmicEmu_etc/32 Models Mira Titan/linear/RUN1/oneh_matterpower.dat")
plot(log10(sim$V1), scrP(sim$V2,sim$V1), type="l",
     xlab=expression(log[10](k)), ylab="script P(k)", ylim=c(0,1.3), xlim = c(-2.1,0),
     main = "LT and simulated spectra")
for (run in 2:10) {
  sim <- read.table(paste0("../CosmicEmu_etc/32 Models Mira Titan/linear/RUN",
                           run,"/oneh_matterpower.dat"))
  lines(log10(sim$V1), scrP(sim$V2,sim$V1), col=run)
}

# From one of the "pow*" folders
for (run in 11:20) {
  if(run == 27 | run == 28) next
  if(run < 20)  folder <- ""
  if(run >= 20) folder <- "_2"
  if(run >= 30) folder <- "_3"
  if(run >= 40) folder <- "_4"
  sim2 <- read.table(paste0("../CosmicEmu_etc/32 Models Mira Titan/pow",folder,
                            ".ic/M0", run, 
                            "/L1300/PM000/analysis/Pow/m0", run, ".pk.ini"))
  lines(log10(sim2$V1)[sim2$V1<=1], 
        scrP(sim2$V2/(Dp[run-10])^2, sim2$V1)[sim2$V1<=1], col=run-10, lwd=3, lty=2)
}

legend(x = "bottom", legend = c("linear theory", "sims"),
      lty = c(1,2), lwd=c(1,3))
