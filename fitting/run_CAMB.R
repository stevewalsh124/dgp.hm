
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

PDF = T
if(PDF) pdf("results/run_CAMB.pdf")
coverages = c()
cov_mat = list()

# Read command line arguments -------------------------------------------------

for(model in 1:32){ # integer 1-32
  print(model)
  deep <- 1 
  
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
  
  # read in and plot all n_lo=15 LR runs (adjusted by h and h^3)
  n_lo = 15
  lr = read.table(paste0("../CosmicEmu_etc/pows/M0",if(model<10){"0"},model,
                         "/L1300/PM001/output/m0",
                         if(model<10){"0"},model,".pk.ini"))
  x_lo = log10(h*lr[which(h*lr[,1]<=x_ub),1])
  # plot(x_lo, scrP(lr[which(h*lr[,1]<=x_ub),2]/h^3, 10^x_lo),
  #      col="gray", type="l", main = model,
  #      ylab = "script P (adj by h,g)",
  #      xlab = "log10(k) (adj by h)")
  y_lo = matrix(NA, length(x_lo), n_lo)
  y_lo[,1] = scrP(lr[which(h*lr[,1]<=x_ub),2]/h^3, 10^x_lo)
  
  for (i in 2:n_lo) {
    if(model==8 & i==9) next
    lr = read.table(paste0("../CosmicEmu_etc/pows/M0",if(model<10){"0"},model,
                           "/L1300/PM0",if(i<10){"0"},i,"/output/m0",
                           if(model<10){"0"},model,".pk.ini"))
    y_lo[,i] = scrP(lr[which(h*lr[,1]<=x_ub),2]/h^3, 10^x_lo)
    # lines(x_lo, y_lo[,i], col="gray")
  }
  if(model==8){
    y_lo = y_lo[,-9]
    n_lo = 14
  }
  
  
  # read in and plot the HR run (adjusted by h and h^3)
  hr = read.table(paste0("../CosmicEmu_etc/pows/M0",if(model<10){"0"},model,
                         "/L2100/HACC000/output/m0",if(model<10){"0"},model,".pk.ini"))
  x_hi = log10(h*hr[which(h*hr[,1]<=x_ub),1])
  y_hi = scrP(hr[which(h*hr[,1]<=x_ub),2]/h^3, 10^x_hi)
  # lines(x_hi, y_hi, col="red")#, lwd=2)
  
  # plot the CAMB data (adjusted by growth^2, h and h^3)
  camb = read.table(paste0("../CosmicEmu_etc/pow64/RUN",model,
                           "/oneh_matterpower.dat"))
  x_camb = log10(camb[,1]*h)
  y_camb = scrP(camb[,2]/h^3*g^2, 10^x_camb)
  # lines(x_camb, y_camb, col="blue1",type="l", lwd=3, lty=3)
  # legend("bottomright", c("camb","lr","hr"), lty=1, col=c("blue","grey","red"))
  
  # Interpolate responses -------------------------------------------------------
  y_hii = approx(x_hi, y_hi, x_lo)$y
  y_cambi = approx(x_camb, y_camb, x_lo)$y
  # plot(x_lo, y_lo[,1], type="l", col="gray", main = "after interpolation")
  # for (i in 2:n_lo) lines(x_lo, y_lo[,i], col="gray")
  # lines(x_lo, y_hii, col="red")
  # lines(x_lo, y_cambi, col="blue")
  
  # remove final value to avoid NA in y_hii
  if(any(is.na(y_cambi) | is.na(y_hii))){
    x_lo = x_lo[-length(x_lo)]
    y_lo = y_lo[-length(x_lo),]
    y_hii = y_hii[-length(x_lo)]
    y_cambi = y_cambi[-length(x_lo)]
  }

  # Scale inputs ----------------------------------------------------------------
  
  # log10 of wavenumber (k) is x
  x <- x_lo
  x <- (x - min(x)) / (max(x) - min(x))
  dx <- deepgp::sq_dist(x)
  
  # Scale responses -------------------------------------------------------------
  # Scale all of these to the emulation space (scrP function from package)
  
  # Get indices for where each data product is deemed unbiased
  lo_ind <- which(x_lo < -0.5)
  hi_ind <- 1:length(x_lo)
  pt_ind <- 1:length(x_lo)
  hi_only <- hi_ind[which(!(hi_ind %in% lo_ind))]
  
  # Adjust the low-res precision info based on the scaling of the response
  var_lo <- apply(y_lo, 1, var)
  # plot(x_lo, log(var_lo),
  #      main = "log k vs log Var(P(k))")
  lm_var <- lm(log(var_lo) ~ x_lo)# + I(x_lo^2))
  # lines(x_lo, lm_var$fitted.values, col="blue", lwd=3)
  prec = as.numeric(1/exp(lm_var$fitted.values))
  prec_adj = 2#3.725626
  prec_lo = ifelse(1:length(x_lo) %in% lo_ind, prec, 0)
  prec_hi = ifelse(1:length(x_lo) %in% hi_ind, prec, 0) * prec_adj
  
  # Obtain low res average
  y_loa <- rowMeans(y_lo)
  
  # Get a weighted average across low, high and pert theory
  # See "Weighted average from multiple computer experiments"
  # in Walsh dissertation: 3.7.1 Appendix E: Derivations
  y_avg <- (1 / (prec_lo*n_lo+prec_hi)) * (prec_lo*n_lo*y_loa + prec_hi*y_hii)
  
  # Scale the responses
  mean_y <- mean(y_avg)
  sd_y <- sd(y_avg)
  y_avg <- (y_avg - mean_y) / sd_y
  for (i in 1:ncol(y_lo)) y_lo[, i] <- (y_lo[, i] - mean_y) / sd_y 
  y_hii <- (y_hii - mean_y) / sd_y 
  y_cambi <- (y_cambi - mean_y) / sd_y 
  
  prec_lo_sz <- prec_lo * sd_y^2
  sd_lo_sz <- sqrt(1 / prec_lo_sz)
  
  # plot(x_lo, y_lo[,1], type="l", col="gray", main = "after scaling")
  # for(i in 2:n_lo) lines(x_lo, y_lo[,i], type="l", col="gray")
  # lines(x_lo, y_hii, col="red")
  # lines(x_lo, y_cambi, col="blue")
  # lines(x_lo, y_avg, col="green", lty=2)
  
  # Get Sigma_hat ---------------------------------------------------------------
  
  # Get smoothed mean and subtract it from the low res runs
  loess_fit <- loess(y_avg ~ x, span = .5)
  y_lo <- y_lo - loess_fit$fitted
  
  # # Optimize kernel hyperparameters for Matern kernel of low res
  # params <- opt_matern(dx[lo_ind, lo_ind], y_lo[lo_ind, ], sd_lo_sz[lo_ind])
  # print(params)
  # Matern_hat <- deepgp:::Matern(dx[lo_ind, lo_ind], params$tau2_hat,
  #                               params$theta_hat, 1e-8, 2.5)
  # 
  # # Create precision block matrices (blocks correspond to pert, lo, high)
  # # block1 <- diag(rep(10000, length(pt_ind)) * sd_y^2)
  # # Create a precision matrix for the low-res portion only
  # block2lo <- solve(diag(sd_lo_sz[lo_ind]) %*% Matern_hat %*% diag(sd_lo_sz[lo_ind]))
  # # Pad the pert and hi-res portions with zeros
  # block2 <- as.matrix(Matrix::bdiag(#diag(rep(0,length(pt_ind))),
  #                                   block2lo, diag(rep(0, length(hi_only)))))
  # block3 <- diag(prec_adj*prec_lo_sz)
  # Sigma_hat <- solve(block2 + block3)
  Sigma_hat = diag(apply(y_lo, 1, var))/(n_lo+prec_adj)
  
  # Run MCMC --------------------------------------------------------------------
  
  if (deep) {
    # load in initialized estimates for warping and hyperparameters for fit
    # w_0 <- read.csv("results/w0.csv")[[1]]
    # params0 <- read.csv("results/params0.csv")
    fit <- fit_two_layer_hm(x, y_avg, nmcmc = 10000, #w_0 = w_0, 
                            #theta_y_0 = params0$theta_y0,
                            #theta_w_0 = params0$theta_w0,
                            Sigma_hat = Sigma_hat)
  } else {
    fit <- fit_one_layer_hm(x, y_avg, nmcmc = 10000, Sigma_hat = Sigma_hat)
  }
  
  # plot(fit) # optionally investigate trace plots
  fit <- trim(fit, 2000, 4)
  fit <- est_true(fit)
  
  coverages[model] = round(mean(fit$lb < y_cambi & fit$ub > y_cambi), 3)
  print(coverages[model])
  cov_mat[[model]] = (fit$lb < y_cambi & fit$ub > y_cambi)

  plot(x_lo, y_lo[,1], type="l", col="gray90", main = paste(model, "cover",
       coverages[model]))
  abline(h=0, col="black", lty=2)
  for (i in 2:n_lo) lines(x_lo, y_lo[,i], col="gray90")
  lines(x_lo, y_hii - loess_fit$fitted, col="red")
  lines(x_lo, y_cambi - loess_fit$fitted, col="blue")
  lines(x_lo, y_avg - loess_fit$fitted, col="green")
  lines(x_lo, fit$m - loess_fit$fitted, col="orange")
  lines(x_lo, fit$ub - loess_fit$fitted, col="orange")
  lines(x_lo, fit$lb - loess_fit$fitted, col="orange")
  legend("topright", c("camb","lr","hr","avg","dgp"), lty=1, 
         col=c("blue","grey","red","green","orange"))
  
  # Unscale results before storing
  results <- data.frame(x = x_lo, 
                        y = y_avg * sd_y + mean_y, 
                        m = fit$m * sd_y + mean_y, 
                        ub = fit$ub * sd_y + mean_y, 
                        lb = fit$lb * sd_y + mean_y,
                        ubb = fit$ubb * sd_y + mean_y, 
                        lbb = fit$lbb * sd_y + mean_y)
  write.csv(results, paste0("results/", ifelse(deep, "dgp", "gp"), "_",
                            model, "camb.csv"), row.names = FALSE)
  # if(model_name == "M001") save(fit, file = "results/fit_M001camb.rda")
}
if(PDF) dev.off()

mean(coverages)
colMeans(matrix(unlist(cov_mat),32,byrow = T))
