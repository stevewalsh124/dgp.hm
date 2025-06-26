# collect the 32 posterior means from CAMB model fits
post_means <- data.frame(matrix(NA, nrow = 400, ncol = 32))
for (i in 1:32) {
  model_name <- i
  results <- read.csv(paste0("results/CAMB/dgp_", model_name, "postmean_int.csv"))[,2]
  post_means[, i] <- results
  colnames(post_means)[i] <- model_name
}
write.csv(post_means, "results/CAMB/post_means_int_CAMB.csv", row.names = FALSE)

x_unif = read.csv(paste0("results/CAMB/dgp_", model_name, "postmean_int.csv"))[,1]
write.csv(x_unif, "results/CAMB/xcamb_int.csv", row.names = F)


library(dgp.hm)
library(pracma)

# for comparison, collect the 64 inf-res runs for all 64 CAMB model runs
inf_res_ints <- data.frame(matrix(NA, nrow = 400, ncol = 64))
for (model in 1:64) {
  camb = read.table(paste0("../CosmicEmu_etc/pow64/RUN",model,
                           "/oneh_matterpower.dat"))
  h = read.csv("../CosmicEmu_etc/pow64/cambDesigns_32x6x2.csv")[model,2]
  x_camb = log10(camb[,1]*h)
  y_camb = scrP(camb[,2]/h^3, 10^x_camb)
  inf_res_ints[,model] <- interp1(x_camb, y_camb, x_unif, method = "cubic")
}
write.table(inf_res_ints, file = "results/CAMB/inf_res_int.txt",
            row.names = FALSE, col.names = FALSE)

