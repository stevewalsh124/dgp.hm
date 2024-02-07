
library(dgp.hm)

# Get inputs
load("../Mira-Titan-IV-Data/precision_and_indexes.Rdata") # loads k
x <- log10(k)
x <- (x - min(x)) / (max(x) - min(x))

# Get outputs
model <- 1
model_name <- paste0(if(model < 100) {"0"}, if (model < 10) {"0"}, model)
pk2 <- read.table(paste0("../Mira-Titan-IV-data/Mira-Titan-2021/STEP499/pk_M",
                         model_name,"_test.dat"))
n <- nrow(pk2)
precs_lo <- ifelse(1:n %in% index_list$lowres.ix, prec_lowres, 0) * 16
precs_hi <- ifelse(1:n %in% index_list$highres.ix, prec_highres, 0)
precs_pt <- ifelse(1:n %in% index_list$pert.ix, 10000, 0)
Lam_lo <- diag(precs_lo)
Lam_hi <- diag(precs_hi)
Lam_pt <- diag(precs_pt)
Lam_z <- Lam_pt + Lam_lo + Lam_hi
precs <- diag(Lam_z)
Lam_zi <- solve(Lam_z)
Y <- t(apply(t(pk2[, 3:18]), 1, function(x) log10(x*(k^1.5)/(2*pi^2))))
Y <- ifelse(is.infinite(Y), 0, Y)
y_lo_avg <- colMeans(Y)
y_pt <- log10(pk2[, 2]*k^1.5/(2*pi^2))
y_pt <- ifelse(is.infinite(y_pt), 0, y_pt)
y_hi <- log10(pk2[, 19]*(k^1.5)/(2*pi^2))
y_hi <- ifelse(is.infinite(y_hi), 0, y_hi)
mu_z <- Lam_zi %*% (Lam_pt %*% y_pt + Lam_lo %*% t(t(y_lo_avg)) + Lam_hi %*% y_hi)
y_avg <- c((mu_z - mean(mu_z)) / sd(mu_z))

# Get Sigma_hat
Sigma_hat <- cov(Y) / 16

fitcov <- fit_one_layer_SW(x, y_avg, nmcmc = 5000, true_g = 1e-6, Sigma_hat = Sigma_hat)
  # could fix true_g = 1e-6
plot(fitcov) # from deepgp package


