
###############################################################################
# For a particular model, plot the Sigma_hat matrix.
#
# Specifications:
#    model: model number (0-116)
#
###############################################################################

library(dgp.hm)
library(zoo) 
library(Matrix) 
library(fields)

# Get model data --------------------------------------------------------------

model <- 1 # integer 0-116 (0, 112:116 are testing)
if (model <= 111) {
  model_name <- paste0("M", if (model < 100) {"0"}, if (model < 10) {"0"}, model) 
} else {
  test_names <- c("E001", "E002", "E003", "E009", "E010")
  model_name <- test_names[model - 111]
}

file_name <- paste0("../Mira-Titan-IV-data/Mira-Titan-2021/STEP499/pk_", 
                    model_name, "_test.dat")
pk2 <- read.table(file_name)
n <- nrow(pk2)
nrun <- 16 # number of low res runs

# Get precisions --------------------------------------------------------------

# Load precision info (k, prec_highres, prec_lowres, index_list)
load("../Mira-Titan-IV-Data/precision_and_indexes.Rdata")

# Get precision info for each data type (low-res, hi-res, pert theory)
precs_lo <- ifelse(1:n %in% index_list$lowres.ix, prec_lowres, 0) * nrun
precs_hi <- ifelse(1:n %in% index_list$highres.ix, prec_highres, 0)
precs_pt <- ifelse(1:n %in% index_list$pert.ix, 10000, 0)

# Smooth the precision information (to remove steps)
Lam_lo <- rollmean(precs_lo, k = 10, fill = "extend")
Lam_hi <- rollmean(precs_hi, k = 10, fill = "extend")
Lam_pt <- rollmean(precs_pt, k = 10, fill = "extend")

# Precision matrix for weighted average, mu_z
Lam_z <- Lam_pt + Lam_lo + Lam_hi

# Get inputs ------------------------------------------------------------------

# log10 of wavenumber (k) is x
x <- log10(k)
x <- (x - min(x)) / (max(x) - min(x))
dx <- deepgp::sq_dist(x)

# Get response ----------------------------------------------------------------

y_pt <- pk2[, 2]
y_lo <- as.matrix(pk2[, 3:18]) # each column is one run
y_hi <- pk2[, 19]

# Scale all of these to the emulation space (scrP function from package)
y_pt <- scrP(y_pt, k)
y_lo <- scrP(y_lo, k)
y_hi <- scrP(y_hi, k)

# Obtain low res average
y_lra <- rowMeans(y_lo)

# Get a weighted average across low, high and pert theory
# See "Weighted average from multiple computer experiments"
# in Walsh dissertation: 3.7.1 Appendix E: Derivations
y_avg <- (1 / Lam_z) * (Lam_pt * y_pt + Lam_lo * y_lra + Lam_hi * y_hi)

# Adjust effective sample size accordingly
# 16 low-res, and hi-res is 3.72x more precise than low-res
hi_wt <- unique(prec_highres / prec_lowres)[1]
nrun <- nrun + hi_wt

# Scale the responses
mean_y <- mean(y_avg)
sd_y <- sd(y_avg)
y_avg <- (y_avg - mean_y) / sd_y
for (i in 1:ncol(y_lo)) y_lo[, i] <- (y_lo[, i] - mean_y) / sd_y 

# Get Sigma_hat ---------------------------------------------------------------

# Get indices
lo_ind <- index_list$lowres.ix
hi_ind <- index_list$highres.ix
pt_ind <- index_list$pert.ix
hi_only <- hi_ind[which(!(hi_ind %in% lo_ind))]

# Adjust the precision info based on the scaling of the response
prec <- Lam_z * sd_y^2
sdd <- sqrt(1 / prec)

# Get smoothed mean and subtract it from the low res runs
loess_fit <- loess(y_avg ~ x, span = 0.15)
y_lo <- y_lo - loess_fit$fitted

# OLD VERSION 
params <- opt_matern(dx[lo_ind, lo_ind], y_lo[lo_ind, ], sdd[lo_ind])
params
Matern_hat <- params$tau2_hat * (geoR::matern(sqrt(dx[lo_ind, lo_ind]), 
                                              phi = params$phi_hat, 
                                              kappa = params$kappa_hat) + 
                                   diag(params$g_hat, length(lo_ind)))

# NEW VERSION
params <- opt_matern2(dx[lo_ind, lo_ind], y_lo[lo_ind, ], sdd[lo_ind])
params
Matern_hat2 <- deepgp:::Matern(dx[lo_ind, lo_ind], params$tau2_hat, 
                              params$theta_hat, 1e-8, 1.5)

par(mfrow = c(1, 2))
image.plot(Matern_hat, main = "Old version, est nugget and kappa")
image.plot(Matern_hat2, main = "New version, g = 1e-8, v = 1.5")

# Create block matrix (blocks correspond to pert, lo, high)
block1 <- (1 / sdd[pt_ind]) * (1/10000)
block2 <-  diag(sdd[lo_ind]) %*% Matern_hat %*% diag(sdd[lo_ind])
block3 <- (1/sdd[hi_only]) * (1 / (precs_hi[hi_only] * sd_y^2))
Sigma_hat <- as.matrix(Matrix::bdiag(diag(block1), block2, diag(block3)))

# Image plot ------------------------------------------------------------------

image.plot(Sigma_hat / nrun)

