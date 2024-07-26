
# Compare my timings to Steve's
annie <- read.csv("results/sims_1_4_5.csv")
steve <- read.csv("results/STEVE_sims_1_4_5.csv")
annie[, 6:9]
steve[, 6:9]

r <- 5 
par(mfrow = c(2, 4))
for (func in 1:2) {
  for (setting in 1:4) {
    filename <- paste0("results/sims_", func, "_", setting, "_", r, ".csv")
    results <- read.csv(filename)[, 2:5]
    boxplot(results, ylab = "MSE", 
            main = paste0("Function ", func, " Setting ", setting))
  }
}
