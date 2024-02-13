
###############################################################################
# For a particular model, plot the fitted mean and variance alongside the
# original data.
#
# Specifications:
#    model: model number (0-116)
#
###############################################################################

# Get model data
model <- 1 
if (model <= 111) {
  model_name <- paste0("M", if(model < 100) {"0"}, if (model < 10) {"0"}, model) 
} else {
  test_names <- c("E001", "E002", "E003", "E009", "E010")
  model_name <- test_names[model - 111]
}

# Get original response data
load("../Mira-Titan-IV-Data/precision_and_indexes.Rdata")
x <- log10(k)
precs_lo <- ifelse(1:n %in% index_list$lowres.ix, prec_lowres, 0) * nrun
precs_hi <- ifelse(1:n %in% index_list$highres.ix, prec_highres, 0)
precs_pt <- ifelse(1:n %in% index_list$pert.ix, 10000, 0)
Lam_lo <- rollmean(precs_lo, k = 10, fill = "extend")
Lam_hi <- rollmean(precs_hi, k = 10, fill = "extend")
Lam_pt <- rollmean(precs_pt, k = 10, fill = "extend")
Lam_z <- Lam_pt + Lam_lo + Lam_hi
y_pt <- pk2[, 2]
y_lo <- as.matrix(pk2[, 3:18]) # each column is one run
y_hi <- pk2[, 19]
y_pt <- scrP(y_pt, k)
y_lo <- scrP(y_lo, k)
y_hi <- scrP(y_hi, k)
y_lra <- rowMeans(y_lo)
y_avg <- (1 / Lam_z) * (Lam_pt * y_pt + Lam_lo * y_lra + Lam_hi * y_hi)

# Get fitted values
fit <- read.csv(paste0("../fitting/results/dgp_", model_name, ".csv"))

# Plot fitted values over the original data
plot(fit$x, fit$y, type = "l", lwd = 1.5, lty = 2, 
     ylim = range(c(fit$y, fit$lbb, fit$ubb)))
for (i in 1:16) lines(fit$x, y_lo[, i], col = "gray")
lines(fit$x, y_hi, lwd = 1.5, col = "red")
lines(fit$x, fit$m, col = "blue")
lines(fit$x, fit$lb, col = "blue", lty = 2)
lines(fit$x, fit$ub, col = "blue", lty = 2)

# TODO: add cosmic emu?


