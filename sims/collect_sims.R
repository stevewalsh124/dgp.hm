###############################################################################
# This script compares results for the "bake-off" of various surrogate models 
# under 2 different test functions (f1 and f2) and 4 different noise settings 
# for a specified random seed. All models are tasked with estimating the noise.
###############################################################################

r <- 5 
par(mfrow = c(1,2))
for (func in 1:2) {
  for (setting in 1:4) {
    filename <- paste0("results/sims_", func, "_", setting, "_", r, ".csv")
    # Compare MSEs
    results <- read.csv(filename)[, c(3,5,4)]
    boxplot(results, ylab = "MSE",
            main = paste0("MSEs: Fn ", func, " Stg ", setting),
            names = c("dgp.hm","hetGP","deepGP"), col = c("green4","gold","gray"))
    # Compare log scores
    results <- read.csv(filename)[, c(7,9,8)]
    boxplot(results, ylab = "log score",
            main = paste0("log scores: Fn ", func, " Stg ", setting),
            names = c("dgp.hm","hetGP","deepGP"), col = c("green4","gold","gray"))
    # # Compare timing
    # results <- read.csv(filename)[, c(11,13,12)]
    # boxplot(results, ylab = "Timing",
    #         main = paste0("timing: Fn ", func, " Stg ", setting),
    # names = c("dgp.hm","hetGP","deepGP"), col = c("green4","gold","gray"))
  }
}
