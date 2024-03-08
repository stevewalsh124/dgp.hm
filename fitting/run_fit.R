
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
file_name <- paste0("../Mira-Titan-IV-Data/Mira-Titan-2021/STEP499/pk_", 
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
precs_pt <- ifelse(1:n %in% index_list$pert.ix, 100000, 0)

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

# Adjust the low-res precision info based on the scaling of the response
prec_lo_sz <- Lam_lo * sd_y^2
sd_lo_sz <- sqrt(1 / prec_lo_sz)

# Get smoothed mean and subtract it from the low res runs
loess_fit <- loess(y_avg ~ x, span = 0.15)
y_lo <- y_lo - loess_fit$fitted

# Optimize kernel hyperparameters for Matern kernel of low res 
params <- opt_matern(dx[lo_ind, lo_ind], y_lo[lo_ind, ], sd_lo_sz[lo_ind])
Matern_hat <- deepgp:::Matern(dx[lo_ind, lo_ind], params$tau2_hat, 
                              params$theta_hat, 1e-8, 2.5)

# Create precision block matrices (blocks correspond to pert, lo, high)
block1 <- diag(Lam_pt * sd_y^2)
# Create a precision matrix for the low-res portion only
block2lo <- solve(diag(sd_lo_sz[lo_ind]) %*% Matern_hat %*% diag(sd_lo_sz[lo_ind]))
# Pad the pert and hi-res portions with zeros
block2 <- as.matrix(Matrix::bdiag(diag(rep(0,length(pt_ind))), 
                                  block2lo, diag(rep(0, length(hi_only)))))
block3 <- diag(Lam_hi * sd_y^2)
Sigma_hat <- solve(block1 + block2 + block3)

# Run MCMC --------------------------------------------------------------------

if (deep) {
  # load in initialized estimates for warping and hyperparameters for fit
  w_0 <- read.csv("results/w0.csv")[[1]]
  params0 <- read.csv("results/params0.csv")
  fit <- fit_two_layer_SW(x, y_avg, nmcmc = 5000, w_0 = w_0, 
                          theta_y_0 = params0$theta_y0,
                          theta_w_0 = params0$theta_w0,
                          Sigma_hat = Sigma_hat)
} else {
  fit <- fit_one_layer_SW(x, y_avg, nmcmc = 5000, Sigma_hat = Sigma_hat)
}

# plot(fit) # optionally investigate trace plots
fit <- trim(fit, 2500, 5)
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
