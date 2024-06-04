###############################################################################
# This script fits the Bayesian hierarchical model for simulated data under
# two different functions: f1 and f2
###############################################################################

library(dgp.hm)

r <- 50
j <- 1

results <- list()

for(model in 1:2){
  for(setting in 1:3){
    result <- read.csv(paste0("results/sims_",model,"_",setting,"_",r,".csv"))$x
    results[[j]] <- result
    j <- j + 1
  }
}

boxplot(results, names = c("1,A","1,B","1,C","2,A","2,B","2,C"), 
        ylab="Average MSE", main = "DGP-HM MSE across models/settings")

# Boxplot for f1
boxplot(results[1:3], ylim=c(0,.018), main="mses: function 1",
        names = c("1,A","1,B","1,C"), ylab="MSE")
grid(nx = NULL, ny = NULL,
     col = "#ebebeb", lwd = 2, lty=1)
boxplot(results[1:3], add = TRUE, names=F)

# Boxplot for f2
boxplot(results[4:6], ylim=c(0,.005), main="mses: function 2",
        names = c("2,A","2,B","2,C"), ylab="MSE")
grid(nx = NULL, ny = NULL,
     col = "#ebebeb", lwd = 2, lty=1)
boxplot(results[4:6], add = TRUE, names = F)
