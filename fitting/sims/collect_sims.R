
func <- 1
setting <- 1
r <- 5 

par(mfrow = c(1, 2))
for (setting in 1:2) {
  filename <- paste0("results/sims_", func, "_", setting, "_", r, ".csv")
  results <- read.csv(filename)[, -1]
  boxplot(results)
}
