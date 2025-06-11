
###############################################################################
# This script fits the Bayesian hierarchical model for a particular CAMB model
# and writes the results to a csv file.
#
# Command line arguments:
#    model: number of particular model (1-32)
#    deep: indicator for whether to fit a DGP or GP (1 = deep, 0 = not deep)
# 
# Outputs:
#    writes csv file of predicted values to the results folder
#
###############################################################################

library(dgp.hm)
library(zoo) # rollmean
library(Matrix) # bdiag
library(pracma) #interp1

# keep track of coverages for CAMB based on DGP fit
coverages = c()
cov_mat = list()

for(model in 1:32){ # integer 1-32
  print(model)
  deep = 1 
  
  # Read command line arguments -------------------------------------------------
  
  args <- commandArgs(TRUE)
  if(length(args) > 0) 
    for(i in 1:length(args)) 
      eval(parse(text = args[[i]]))
  
  cat("model is ", model, "\n")
  if (deep) cat("model is deep \n") else cat("model is NOT deep \n")
  
  # Load data for a particular model --------------------------------------------
  # CAMB has runs 1-64; L1300 & L2100 has runs 1-32
  
  # growth values; need to multiply CAMB by growth^2 (only 32 tho...)
  g = read.table("../CosmicEmu_etc/pows/growth.dat")[model,2]
  # hubble (h) values; mult k by h, divide P(k) by h^3
  h = read.csv("../CosmicEmu_etc/pow64/cambDesigns_32x6x2.csv")[model,2]
  
  # upper bound to cut off x values (orig scale, not log10)
  x_ub = 0.2
  
  # read in all n_lo=15 low-res (L1300) runs (adjusted by h and h^3)
  n_lo = 15
  lr = read.table(paste0("../CosmicEmu_etc/pows/M0",if(model<10){"0"},model,
                         "/L1300/PM001/output/m0",
                         if(model<10){"0"},model,".pk.ini"))
  x_lo = log10(h*lr[which(h*lr[,1]<=x_ub),1])
  y_lo = matrix(NA, length(x_lo), n_lo)
  y_lo[,1] = scrP(lr[which(h*lr[,1]<=x_ub),2]/h^3/g^2, 10^x_lo)
  
  for (i in 2:n_lo) {
    if(model==8 & i==9) next
    lr = read.table(paste0("../CosmicEmu_etc/pows/M0",if(model<10){"0"},model,
                           "/L1300/PM0",if(i<10){"0"},i,"/output/m0",
                           if(model<10){"0"},model,".pk.ini"))
    y_lo[,i] = scrP(lr[which(h*lr[,1]<=x_ub),2]/h^3/g^2, 10^x_lo)
  }
  
  # only model 8 is missing the 9th run
  if(model==8){ 
    y_lo = y_lo[,-9]
    n_lo = 14
  }
  
  # read in the high-res (L2100) run (adjusted by h and h^3)
  hr = read.table(paste0("../CosmicEmu_etc/pows/M0",if(model<10){"0"},model,
                         "/L2100/HACC000/output/m0",if(model<10){"0"},model,".pk.ini"))
  x_hi = log10(h*hr[which(h*hr[,1]<=x_ub),1])
  y_hi = scrP(hr[which(h*hr[,1]<=x_ub),2]/h^3/g^2, 10^x_hi)

  # read in the CAMB data (adjusted by growth^2, h and h^3)
  camb = read.table(paste0("../CosmicEmu_etc/pow64/RUN",model,
                           "/oneh_matterpower.dat"))
  x_camb = log10(camb[,1]*h)
  y_camb = scrP(camb[,2]/h^3, 10^x_camb)

  # Interpolate responses -------------------------------------------------------
  # Create x's of [-2.5, -2.2] for lower-res version of x_camb, use with x_lo
  x = c(seq(-2.5, round(min(x_lo), 1), by = 0.1), x_lo)
  n_camb_only = length(seq(-2.5, round(min(x_lo), 1), by = 0.1))
  y_hii = approx(x_hi, y_hi, x)$y
  y_cambi = approx(x_camb, y_camb, x)$y
  y_loi = matrix(NA, length(x), n_lo)
  for (i in 1:n_lo) y_loi[,i] = approx(x_lo, y_lo[,i], x)$y

  # if NAs, change to 0s 
  y_hii = ifelse(is.na(y_hii), 0, y_hii)
  y_loi = ifelse(is.na(y_loi), 0, y_loi)
  y_cambi = ifelse(is.na(y_cambi), 0, y_cambi)

  # Scale based on precision indices and values ---------------------------------
  
  # Adjust the low-res precision info based on the scaling of the response
  var_lo <- apply(y_lo, 1, var)
  lm_var <- lm(log(var_lo) ~ x_lo)# + I(x_lo^2))

  # Combine 1e8 prec (for camb) with log-log precs for lo and hi
  prec = c(rep(1e8, n_camb_only), as.numeric(1/exp(lm_var$fitted.values)))
  prec_adj = 3.725626
  
  # Get indices for where each data product is deemed unbiased
  lo_ind <- which(x < log10(x_ub) & x >= min(x_lo))
  hi_ind <- which(x < log10(x_ub) & x >= min(x_lo))
  camb_ind <- which(x <= min(x_lo))
  hi_only <- hi_ind[which(!(hi_ind %in% lo_ind))]
  
  # Smooth the precision information (to remove steps)
  prec_lo = pmax(0, smooth.spline(rollmean(ifelse(1:length(x) %in% lo_ind, prec, 0),
                                           k = 3, fill = "extend"))$y)
  prec_hi = pmax(0, smooth.spline(rollmean(ifelse(1:length(x) %in% hi_ind, prec, 0) * prec_adj,
                                           k = 3, fill = "extend"))$y)
  prec_camb = pmax(0, smooth.spline(rollmean(ifelse(1:length(x) %in% camb_ind, 1e8, 0),
                                             k = 3, fill = "extend"))$y)
  prec_avg = (prec_lo*n_lo + prec_hi + prec_camb)
  
  # Scale inputs ----------------------------------------------------------------
  
  # log10 of wavenumber (k) is x
  # use camb for [-4.2,-2.1], and others when > -2.1
  b <- max(x)
  a <- min(x)
  x <- (x - a) / (b - a)
  dx <- deepgp::sq_dist(x)
  
  # Scale responses -------------------------------------------------------------
  # Scale all of these to the emulation space (scrP function from package)
  
  # Obtain low res average
  y_loa <- rowMeans(y_loi)
  
  # Get a weighted average across low, high and camb theory
  # See "Weighted average from multiple computer experiments"
  # in Walsh dissertation: 3.7.1 Appendix E: Derivations
  y_avg <- (1 / prec_avg) * 
    (prec_camb*y_cambi + prec_lo*n_lo*y_loa + prec_hi*y_hii)
  
  # Scale the responses
  mean_y <- mean(y_avg)
  sd_y <- sd(y_avg)
  y_avg <- (y_avg - mean_y) / sd_y
  for (i in 1:ncol(y_loi)) y_loi[, i] <- (y_loi[, i] - mean_y) / sd_y 
  y_hii <- (y_hii - mean_y) / sd_y 
  y_cambi <- (y_cambi - mean_y) / sd_y 
  
  prec_avg_sz <- prec_avg * sd_y^2
  
  # Get Sigma_hat ---------------------------------------------------------------
  
  # Get smoothed mean and subtract it from the low res runs
  loess_fit <- loess(y_avg ~ x, span = .5)
  # y_lo <- y_lo - loess_fit$fitted
  # lines(x, loess_fit$fitted, col="green", lty=3)
  
  # # Optimize kernel hyperparameters for Matern kernel of low res
  # prec_lo_sz <- prec_lo * sd_y^2
  # sd_lo_sz <- sqrt(1 / prec_lo_sz)
  # params <- opt_matern(dx[lo_ind, lo_ind], y_lo[lo_ind, ], sd_lo_sz[lo_ind])
  # print(params)
  # Matern_hat <- deepgp:::Matern(dx[lo_ind, lo_ind], params$tau2_hat,
  #                               params$theta_hat, 1e-8, 2.5)
  # 
  # # Create precision block matrices (blocks correspond to camb, lo, high)
  # # block1 <- diag(rep(1e8, length(camb_ind)) * sd_y^2)
  # # Create a precision matrix for the low-res portion only
  # block2lo <- solve(diag(sd_lo_sz[lo_ind]) %*% Matern_hat %*% diag(sd_lo_sz[lo_ind]))
  # # Pad the camb and hi-res portions with zeros
  # block2 <- as.matrix(Matrix::bdiag(#diag(rep(0,length(camb_ind))),
  #                                   block2lo, diag(rep(0, length(hi_only)))))
  # block3 <- diag(prec_adj*prec_lo_sz)
  # Sigma_hat <- solve(block2 + block3)
  Sigma_hat = diag(1/prec_avg_sz)
  
  # Run MCMC --------------------------------------------------------------------
  
  if (deep) {
    # load in initialized estimates for warping and hyperparameters for fit
    # w_0 <- read.csv("results/w0.csv")[[1]]
    # params0 <- read.csv("results/params0.csv")
    fit <- fit_two_layer_hm(x, y_avg, nmcmc = 15000, #w_0 = w_0, 
                            # theta_y_0 = params0$theta_y0,
                            # theta_w_0 = params0$theta_w0,
                            # settings = list(pmx=T),
                            Sigma_hat = Sigma_hat)
  } else {
    fit <- fit_one_layer_hm(x, y_avg, nmcmc = 10000, Sigma_hat = Sigma_hat)
  }
  
  # plot(fit) # optionally investigate trace plots
  fit <- trim(fit, 10000, 4)
  fit <- est_true(fit)
  
  coverages[model] = round(mean(fit$lb < y_cambi & fit$ub > y_cambi), 3)
  print(coverages[model])
  cov_mat[[model]] = (fit$lb < y_cambi & fit$ub > y_cambi)
  
  # interpolate posterior mean to a high-resolution uniform grid
  x_unif = seq(0, 1, length.out = 400)
  # # plot with loess avg subtracted...
  # fitm_int1 <- interp1(x, fit$m - loess_fit$fitted, x_unif, method = "cubic")
  # lines(x_unif, fitm_int1, col="purple",lwd=2, lty=3)
  # ...but save a version without loess subtracted
  fitm_int1 <- interp1(x, fit$m, x_unif, method = "cubic")
  
  # Unscale results before storing
  results <- data.frame(x = x*(b-a)+a, 
                        y = y_avg * sd_y + mean_y, 
                        m = fit$m * sd_y + mean_y, 
                        ub = fit$ub * sd_y + mean_y, 
                        lb = fit$lb * sd_y + mean_y,
                        ubb = fit$ubb * sd_y + mean_y, 
                        lbb = fit$lbb * sd_y + mean_y)
  results_int = data.frame(x = x_unif*(b-a)+a, y = fitm_int1*sd_y+mean_y)
  write.csv(results, paste0("results/CAMB/", ifelse(deep, "dgp", "gp"), "_",
                            model, "camb.csv"), row.names = FALSE)
  write.csv(results_int, paste0("results/CAMB/", ifelse(deep, "dgp", "gp"), "_",
                                model, "postmean_int.csv"), row.names = FALSE)
  if(model == 1) save(fit, file = "results/CAMB/fit_1camb.rda")
  if(model == 1) save.image(file = "results/CAMB/fit_1camb_image.rda")
}

mean(coverages)
# colMeans(matrix(unlist(cov_mat), 4, byrow = T))
