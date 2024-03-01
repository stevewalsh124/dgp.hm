
###############################################################################
# This script grabs the posterior mean from each fit.
#
# Models 1-111 are saved in "post_means_train.csv".
# Models 0, 112-116 are saved in "post_means_test.csv".
#
###############################################################################

# Training models
post_means <- data.frame(matrix(NA, nrow = 351, ncol = 10)) # CHANGE TO 111 ONCE THEY ARE RUN
for (i in 1:10) {
  model_name <- paste0("M", if (i < 100) {"0"}, if (i < 10) {"0"}, i) 
  results <- read.csv(paste0("results/dgp_", model_name, ".csv"))
  post_means[, i] <- results$m
  colnames(post_means)[i] <- model_name
}
write.csv(post_means, "results/post_means_train.csv", row.names = FALSE)

# Testing models (STILL NEED TO BE RUN FIRST)
#post_means <- data.frame(matrix(NA, nrow = 351, ncol = 7))
#test_names <- c("M000", "E001", "E002", "E003", "E009", "E010")
#for (i in 1:7) {
#  model_name <- test_names[i]
#  results <- read.csv(paste0("results/dgp_", model_name, ".csv"))
#  post_means[, i] <- results$m
#  colnames(post_means)[i] <- model_name
#}
#write.csv(post_means, "post_means_test.csv", row.names = FALSE)
