# Power Spectra Simulations

library(dgp.hm)

# From the "linear" folder
sim <- read.table("../CosmicEmu_etc/32 Models Mira Titan/linear/RUN1/oneh_matterpower.dat")
plot(log10(sim$V1), scrP(sim$V2,sim$V1), type="l",
     xlab=expression(log[10](k)), ylab="script P(k)", ylim=c(-4.7,1.1),
     main = "all 'linear' runs")
for (run in 2:32) {
  sim <- read.table(paste0("../CosmicEmu_etc/32 Models Mira Titan/linear/RUN",
                           run,"/oneh_matterpower.dat"))
  lines(log10(sim$V1), scrP(sim$V2,sim$V1), col=run)
}

# From one of the "pow*" folders
sim2 <- read.table("../CosmicEmu_etc/32 Models Mira Titan/pow.ic/M011/L1300/PM000/analysis/Pow/m011.pk.ini")
head(sim2)
pairs(sim2)
plot(log10(sim2$V1), log10(sim2$V2), main = "orig, no cutoff")
plot(log10(sim2$V1), scrP(sim2$V2, sim2$V1), main = "emu, no cutoff")
plot(log10(sim2$V1)[sim2$V1<=1], scrP(sim2$V2, sim2$V1)[sim2$V1<=1], 
     main = "emu, cutoff at k=1", type="l", xlim=c(-2.25,0),ylim=c(-4.3,-2.6))
for (run in 12:42) {
  if(run == 27 | run == 28) next
  if(run < 20)  folder <- ""
  if(run >= 20) folder <- "_2"
  if(run >= 30) folder <- "_3"
  if(run >= 40) folder <- "_4"
  sim2 <- read.table(paste0("../CosmicEmu_etc/32 Models Mira Titan/pow",folder,
                            ".ic/M0", run, 
                            "/L1300/PM000/analysis/Pow/m0", run, ".pk.ini"))
  lines(log10(sim2$V1)[sim2$V1<=1], scrP(sim2$V2, sim2$V1)[sim2$V1<=1], 
       main = "emu, cutoff at k=1", col=run)
}
