
###############################################################################
# For a particular model, plot the different data types in the emulation space.
#
# Specifications:
#    model: model number (0-116)
#
###############################################################################

library(dgp.hm)
library(zoo)

# Get model data
model <- 1 
if (model <= 111) {
  model_name <- paste0("M", if(model < 100) {"0"}, if (model < 10) {"0"}, model) 
} else {
  test_names <- c("E001", "E002", "E003", "E009", "E010")
  model_name <- test_names[model - 111]
}
file_name <- paste0("../Mira-Titan-IV-data/Mira-Titan-2021/STEP499/pk_", 
                    model_name, "_test.dat")
pk2 <- read.table(file_name)
n <- nrow(pk2)
nrun <- 16

# Get inputs
load("../Mira-Titan-IV-Data/precision_and_indexes.Rdata")
x <- log10(k)
x <- (x - min(x)) / (max(x) - min(x))

# Get precision info 
precs_lo <- ifelse(1:n %in% index_list$lowres.ix, prec_lowres, 0) * nrun
precs_hi <- ifelse(1:n %in% index_list$highres.ix, prec_highres, 0)
precs_pt <- ifelse(1:n %in% index_list$pert.ix, 10000, 0)
Lam_lo <- rollmean(precs_lo, k = 10, fill = "extend")
Lam_hi <- rollmean(precs_hi, k = 10, fill = "extend")
Lam_pt <- rollmean(precs_pt, k = 10, fill = "extend")
Lam_z <- Lam_pt + Lam_lo + Lam_hi

# Get responses
y_pt <- pk2[, 2]
y_lo <- as.matrix(pk2[, 3:18]) # each column is one run
y_hi <- pk2[, 19]
y_pt <- scrP(y_pt, k)
y_lo <- scrP(y_lo, k)
y_hi <- scrP(y_hi, k)
y_lra <- rowMeans(y_lo)
y_avg <- (1 / Lam_z) * (Lam_pt * y_pt + Lam_lo * y_lra + Lam_hi * y_hi)

# Plot data types
matplot(x, y_lo, col = "gray", type = "l", lty = 1, ylim = range(y_pt, y_lo, y_hi),
        xlab = "x", ylab = "y")
lines(x, y_pt, col = "green")
lines(x, y_lra, col = "blue")
lines(x, y_hi, col = "red")
lines(x[index_list$pert.ix], rep(-1, length(index_list$pert.ix)), 
      col = "green", lty = 3)
lines(x[index_list$lowres.ix], rep(-1.1, length(index_list$lowres.ix)), 
      col = "blue", lty = 3)
lines(x[index_list$highres.ix], rep(-1.2, length(index_list$highres.ix)), 
      col = "red", lty = 3)
lines(x, y_avg, col = "orange", lty = 2, lwd = 4)


