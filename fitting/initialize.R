
###############################################################################
# This script fits the deep Bayesian hierarchical model for M001 for many 
# MCMC iterations and saves the resulting hyperparameters (w_0, theta_w_0, and
# theta_y_0) to csv files.  They will be used to initialize fits for other 
# models.
# 
# Outputs:
#    writes csv file of final mcmc iteration for w and hyperparameters
#    (theta_y, theta_w)
#
###############################################################################

library(dgp.hm)
library(zoo) # rollmean
library(Matrix) # bdiag

# Load data for M001 ----------------------------------------------------------

model_name <- "M001"

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
params <- opt_matern(dx[lo_ind, lo_ind], y_lo[lo_ind, ], sd_lo_sz[lo_ind],
                     n_multi = 10)
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

fit <- fit_two_layer_hm(x, y_avg, nmcmc = 50000, Sigma_hat = Sigma_hat)

# plot(fit) # optionally investigate trace plots
fit <- trim(fit, 49000, 10)

# Unscale results before storing
w0 <- c(fit$w[, fit$nmcmc])
params0 <- data.frame(theta_y0 = fit$theta_y[fit$nmcmc],
                      theta_w0 = fit$theta_w[fit$nmcmc])

write.csv(params0, paste0("results/params0.csv"), row.names = FALSE)
write.csv(w0, paste0("results/w0.csv"), row.names = FALSE)
