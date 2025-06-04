###############################################################################
# This script compares results for the "bake-off" of various surrogate models 
# under 2 different test functions (f1 and f2) and 4 different noise settings 
# for a specified random seed. All models are tasked with estimating the noise.
###############################################################################

PDF <- F
if(PDF) pdf("collect_sims.pdf")

# color-blind friendly color palette
cbcols <- palette.colors(palette = "Okabe-Ito")
#1 black, 2 orange, 3 skyblue , 4 bluishgreen, 5 yellow 
#6 blue, 7 vermillion,  8 reddishpurple, 9 gray 

###################################################
# Get ylim's for MSE/logScore for both r=5 and 20 #
###################################################

MSE5_range <- MSE20_range <- logS5_range <- logS20_range <- c()

r=5
for (func in 1:2) {
  for (setting in 4:5) {
    filename <- paste0("results/sims_", func, "_", setting, "_", r, ".csv")
    # Compare MSEs
    results <- read.csv(filename)[, c(3,6,5,4)]
    # NOTE: removed one outlier to not skew the ylim
    if(func==1&setting==4) results <- results[-which.max(results[,3]),]
    MSE5_range <- range(MSE5_range, results)
    results <- read.csv(filename)[, c(8,11,10,9)]
    # Note: removes three outliers from setting==5, r==5
    logS5_range <- range(logS5_range, results[results <= 150])
  }
}

r=20
for (func in 1:2) {
  for (setting in 4:5) {
    filename <- paste0("results/sims_", func, "_", setting, "_", r, ".csv")
    # Compare MSEs
    results <- read.csv(filename)[, c(3,6,5,4)]
    MSE20_range <- range(MSE20_range, results)
    results <- read.csv(filename)[, c(8,11,10,9)]
    logS20_range <- range(logS20_range, results)
  }
}

####################
# Boxplots of MSEs #
####################

par(mfrow=c(2,4),oma=c(4,4,1.5,1),mar=c(0,0,0,0))
for (r in c(5,20)) {
  for (setting in 4:5) {
    for (func in 1:2) {
      filename <- paste0("results/sims_", func, "_", setting, "_", r, ".csv")
      # Compare MSEs
      results <- read.csv(filename)[, c(3,6,5,4)]
      boxplot(results, ylab = "MSE", axes=F,
              names = c("dgp.hm","dpc","hetGP","deepGP"), col = cbcols[c(4,5,7,9)],
              ylim = if(r==5){MSE5_range}else{MSE20_range})
      mtext(paste0("Fn ", func, " Stg ", setting, " r=", r), line=-1.5)
      box()
      if(func==1 & setting==4) {axis(2); mtext(side=2,line = 2,"MSE")}
      if(r==20) axis(1, at=1:4, labels = c("dgp.hm","dpc","hetGP","deepGP"))

    }
  }
}

##########################
# Boxplots of log scores #
##########################

par(mfrow=c(2,4),oma=c(4,4,1.5,1),mar=c(0,0,0,0))
for (r in c(5,20)) {
  for (setting in 4:5) {
    for (func in 1:2) {
      filename <- paste0("results/sims_", func, "_", setting, "_", r, ".csv")
      # Compare MSEs
      results <- read.csv(filename)[, c(8,11,10,9)]
      boxplot(results, ylab = "MSE", axes=F,
              names = c("dgp.hm","dpc","hetGP","deepGP"), col = cbcols[c(4,5,7,9)],
              ylim = if(r==5){logS5_range}else{logS20_range})
      mtext(paste0("Fn ", func, " Stg ", setting, " r=", r), line=-1.5)
      box()
      if(func==1 & setting==4) {axis(2); mtext(side=2,line = 2,"log score")}
      if(r==20) axis(1, at=1:4, labels = c("dgp.hm","dpc","hetGP","deepGP"))
    }
  }
}
if(PDF) dev.off()
