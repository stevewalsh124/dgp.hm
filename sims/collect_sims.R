###############################################################################
# This script compares results for the "bake-off" of various surrogate models 
# under 2 different test functions (f1 and f2) and 2 different noise settings 
# for a specified random seed. All models are tasked with estimating the noise.
###############################################################################

JPG <- F

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
    results <- read.csv(filename)[, c(2,5,4,3)]
    # NOTE: removed one outlier to not skew the ylim
    if(func==1&setting==4) results <- results[-which.max(results[,3]),]
    MSE5_range <- range(MSE5_range, results)
    results <- read.csv(filename)[, c(6,9,8,7)]
    # Note: removes three outliers from setting==5, r==5
    logS5_range <- range(logS5_range, results[results <= 150])
  }
}

r=15
for (func in 1:2) {
  for (setting in 4:5) {
    filename <- paste0("results/sims_", func, "_", setting, "_", r, ".csv")
    # Compare MSEs
    results <- read.csv(filename)[, c(2,5,4,3)]
    MSE20_range <- range(MSE20_range, results)
    results <- read.csv(filename)[, c(6,9,8,7)]
    logS20_range <- range(logS20_range, results)
  }
}

####################
# Boxplots of MSEs #
####################

if(JPG) jpeg("../paper/sims_MSE.jpeg", width = 10, height = 4, units = "in", res = 300)

par(mfrow=c(2,4), oma=c(4,4,1.5,1), mar=c(0,0,0,0))

# Settings
functions <- 1:2
settings <- 4:5
r_vals <- c(5, 15)

# Step 1: Precompute y-limits for each r value (i.e., per row)
row_ylims_logs <- list()
for (r in r_vals) {
  all_logs <- NULL
  for (func in functions) {
    for (setting in settings) {
      filename <- paste0("results/sims_", func, "_", setting, "_", r, ".csv")
      results <- read.csv(filename)[, c(2,5,4,3)]
      all_logs <- c(all_logs, unlist(results))
    }
  }
  row_ylims_logs[[as.character(r)]] <- range(all_logs, na.rm = TRUE)
}

# Step 2: Plot 2 rows (r), 4 columns (func/setting combos)
par(mfrow=c(2,4), oma=c(4,4,1.5,1), mar=c(0,0,0,0))

for (r in r_vals) {
  row_ylim <- row_ylims_logs[[as.character(r)]]
  for (func in functions) {
    for (setting in settings) {
      filename <- paste0("results/sims_", func, "_", setting, "_", r, ".csv")
      results <- read.csv(filename)[, c(2,5,4,3)]
      
      boxplot(results, axes = FALSE,
              names = c("DGP.FCO","DPC","hetGP","deepgp"),
              col = cbcols[c(4,5,7,9)],
              ylim = row_ylim)
      
      mtext(paste0("Fn ", func, ", Stg ", LETTERS[setting - 3], ", r=", r), line = -1.5)
      box()
      
      # Axis logic: only show y-axis on leftmost column
      if (func == 1 && setting == 4) {
        axis(2)
        mtext(side = 2, line = 2, "MSE")
      }
      if(r==15) axis(1, at = 1:4, labels = c("DGP.FCO","DPC","hetGP","deepgp"))
    }
  }
}

if(JPG) dev.off()

##########################
# Boxplots of log scores #
##########################

if(JPG) jpeg("../paper/sims_logS.jpeg", width = 10, height = 4, units = "in", res = 300)

par(mfrow=c(2,4), oma=c(4,4,1.5,1), mar=c(0,0,0,0))

# Settings
functions <- 1:2
settings <- 4:5
r_vals <- c(5, 15)

# Step 1: Precompute y-limits for each r value (i.e., per row)
row_ylims_logs <- list()
for (r in r_vals) {
  all_logs <- NULL
  for (func in functions) {
    for (setting in settings) {
      filename <- paste0("results/sims_", func, "_", setting, "_", r, ".csv")
      results <- read.csv(filename)[, c(6,9,8,7)]
      all_logs <- c(all_logs, unlist(results))
    }
  }
  row_ylims_logs[[as.character(r)]] <- range(all_logs, na.rm = TRUE)
}

# Step 2: Plot 2 rows (r), 4 columns (func/setting combos)
par(mfrow=c(2,4), oma=c(4,4,1.5,1), mar=c(0,0,0,0))

for (r in r_vals) {
  row_ylim <- row_ylims_logs[[as.character(r)]]
  for (func in functions) {
    for (setting in settings) {
      filename <- paste0("results/sims_", func, "_", setting, "_", r, ".csv")
      results <- read.csv(filename)[, c(6,9,8,7)]
      
      boxplot(results, axes = FALSE,
              names = c("DGP.FCO","DPC","hetGP","deepgp"),
              col = cbcols[c(4,5,7,9)],
              ylim = row_ylim)
      
      mtext(paste0("Fn ", func, ", Stg ", LETTERS[setting - 3], ", r=", r), line = -1.5)
      box()
      
      # Axis logic: only show y-axis on leftmost column
      if (func == 1 && setting == 4) {
        axis(2)
        mtext(side = 2, line = 2, "log score")
      }
      if(r==15) axis(1, at = 1:4, labels = c("DGP.FCO","DPC","hetGP","deepgp"))
    }
  }
}


if(JPG) dev.off()

