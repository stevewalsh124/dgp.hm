library(dgp.hm)

# CAMB has runs 1-64; L1300 & L2100 has runs 1-32
for (run in 1:32) {

  # growth values; need to multiply CAMB by growth^2 (only 32 tho...)
  g = read.table("../CosmicEmu_etc/pows/growth.dat")[run,2]
  # hubble (h) values; mult k by h, divide P(k) by h^3
  h = read.csv("../CosmicEmu_etc/pow64/cambDesigns_32x6x2.csv")[run,2]
  
  # read in and plot all 15 LR runs (adjusted by h and h^3)
  lr = read.table(paste0("../CosmicEmu_etc/pows/M0",if(run<10){"0"},run,
                         "/L1300/PM001/output/m0",
                         if(run<10){"0"},run,".pk.ini"))
  if(run==1) plot(log10(h*lr[which(h*lr[,1]<=1),1]),
                  scrP(lr[which(h*lr[,1]<=1),2]/h^3,lr[which(h*lr[,1]<=1),1]*h),
                  col="gray", type="l", main = run,
                  ylab = "script P (adj by h,g)",
                  xlab = "log10(k) (adj by h)")# ylim=c(-5,-2))
  if(run!=1) plot(log10(h*lr[which(h*lr[,1]<=1),1]),
                  scrP(lr[which(h*lr[,1]<=1),2]/h^3,lr[which(h*lr[,1]<=1),1]*h),
                  col="gray", main=run, type="l",
                  ylab = "script P (adj by h,g)",
                  xlab = "log10(k) (adj by h)")
  lrs = matrix(NA, 15, 206)
  lrs[1,] = scrP(lr[which(h*lr[,1]<=1),2]/h^3,lr[which(h*lr[,1]<=1),1]*h)
  
  for (i in 2:15) {
    if(run ==8 & i ==9) next
    lr = read.table(paste0("../CosmicEmu_etc/pows/M0",if(run<10){"0"},run,
                           "/L1300/PM0",if(i<10){"0"},i,"/output/m0",
                           if(run<10){"0"},run,".pk.ini"))
    lrs[i,] = scrP(lr[which(h*lr[,1]<=1),2]/h^3,lr[which(h*lr[,1]<=1),1]*h)
    lines(log10(h*lr[which(h*lr[,1]<=1),1]), lrs[i,], col="gray")
  }
  
  # read in and plot the HR run (adjusted by h and h^3)
  hr = read.table(paste0("../CosmicEmu_etc/pows/M0",if(run<10){"0"},run,
                         "/L2100/HACC000/output/m0",if(run<10){"0"},run,".pk.ini"))
  lines(log10(h*hr[which(h*hr[,1]<=1),1]),
        scrP(hr[which(h*hr[,1]<=1),2]/h^3,hr[which(h*hr[,1]<=1),1]*h),
        col="red", lwd=2)
  
  # plot the CAMB data (adjusted by growth^2, h and h^3)
  camb = read.table(paste0("../CosmicEmu_etc/pow64/RUN",run,
                           "/oneh_matterpower.dat"))
  lines(log10(camb[,1]*h), scrP(camb[,2]/h^3*g^2,camb[,1]*h),
        col="blue1",type="l", lwd=3, lty=3)

  legend("bottomright",c("camb","lr","hr"), lty=1,col=c("blue","grey","red"))
}
