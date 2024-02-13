
###############################################################################
# This script fits the Bayesian hierarchical model for a particular cosmology
# and writes the results to a csv file.
#
# Command line arguments:
#    model: number of particular model (0-116)
#    deep: indicator for whether to fit a DGP or GP (1 = deep, 0 = not deep)
# 
# Outputs:
#    writes csv file of predicted values to the results folder
#
###############################################################################

library(dgp.hm)
library(zoo) # rollmean
library(Matrix) # bdiag

# Read command line arguments -------------------------------------------------

model <- 1 # integer 0-116 (0, 112:116 are testing)
deep <- 1 

args <- commandArgs(TRUE)
if(length(args) > 0) 
  for(i in 1:length(args)) 
    eval(parse(text = args[[i]]))

cat("model is ", model, "\n")
if (deep) cat("model is deep \n") else cat("model is NOT deep \n")

# Load data for a particular model --------------------------------------------

if (model <= 111) {
  model_name <- paste0("M", if (model < 100) {"0"}, if (model < 10) {"0"}, model) 
} else {
  test_names <- c("E001", "E002", "E003", "E009", "E010")
  model_name <- test_names[model - 111]
}

# 1st column is k, 2nd is linear pert theory, 3:18 is low-res, 19 is hi-res
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
y_lo <- (y_lo - mean_y) / sd_y 

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

# Optimize kernel hyperparameters for Matern kernel of low res 
params <- opt_matern(dx[lo_ind, lo_ind], y_lo[lo_ind, ], sdd[lo_ind])
Matern_hat <- params$tau2_hat * (geoR::matern(sqrt(dx[lo_ind, lo_ind]), 
                                              phi = params$phi_hat, 
                                              kappa = params$kappa_hat) + 
                                   diag(params$g_hat, length(lo_ind)))

# Create block matrix (blocks correspond to pert, lo, high)
block1 <- (1 / sdd[pt_ind]) * (1/10000)
block2 <-  diag(sdd[lo_ind]^2) * Matern_hat
block3 <- (1/sdd[hi_only]) * (1 / (precs_hi[hi_only] * sd_y^2))
Sigma_hat <- as.matrix(Matrix::bdiag(diag(block1), block2, diag(block3)))

# Run MCMC --------------------------------------------------------------------

if (deep) {
  w_0 <- read.csv("w0_from_mte1_50k.csv")[[1]]
  fit <- fit_two_layer_SW(x, y_avg, nmcmc = 1500, w_0 = w_0,
                          Sigma_hat = Sigma_hat / nrun)
} else {
  fit <- fit_one_layer_SW(x, y_avg, nmcmc = 1500, Sigma_hat = Sigma_hat / nrun)
}

# plot(fit) # optionally investigate trace plots
fit <- trim(fit, 1000, 5)
fit <- est_true(fit)

# Unscale results before storing
results <- data.frame(x = log10(k), 
                      y = y_avg * sd_y + mean_y, 
                      m = fit$m * sd_y + mean_y, 
                      ub = fit$ub * sd_y + mean_y, 
                      lb = fit$lb * sd_y + mean_y,
                      ubb = fit$ubb * sd_y + mean_y, 
                      lbb = fit$lbb * sd_y + mean_y)
write.csv(results, paste0("results/", ifelse(deep, "dgp", "gp"), "_", 
                          model_name, ".csv"), row.names = FALSE)
