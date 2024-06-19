
r <- 5 

par(mfrow = c(2, 4))
for (func in 1:2) {
  for (setting in 1:4) {
    filename <- paste0("results/sims_", func, "_", setting, "_", r, ".csv")
    results <- read.csv(filename)[, -1]
    boxplot(results, ylab = "MSE", 
            main = paste0("Function ", func, " Setting ", setting))
  }
}
