
###############################################################################
# For a particular model, plot the different data types in the emulation space.
# Additionally, plot the fitted mean and variance alongside the original data.
# 
# Specifications:
#    model: model number (0-116)
#
###############################################################################

library(dgp.hm)
library(zoo)

# make JPEGs of plots?
JPG <- FALSE

# Get model data
model <- 1
if (model <= 111) {
  model_name <- paste0("M", if (model < 100) {"0"}, if (model < 10) {"0"}, model) 
} else {
  test_names <- c("E001", "E002", "E003", "E009", "E010")
  model_name <- test_names[model - 111]
}
file_name <- paste0("../Mira-Titan-IV-Data/Mira-Titan-2021/STEP499/pk_", 
                    model_name, "_test.dat")
pk2 <- read.table(file_name)
n <- nrow(pk2)
nrun <- 16

# Get inputs
load("../Mira-Titan-IV-Data/precision_and_indexes.Rdata")
x <- log10(k)

# Get precision info 
precs_lo <- ifelse(1:n %in% index_list$lowres.ix, prec_lowres, 0) * nrun
precs_hi <- ifelse(1:n %in% index_list$highres.ix, prec_highres, 0)
precs_pt <- ifelse(1:n %in% index_list$pert.ix, 10^8, 0)
Lam_lo <- rollmean(precs_lo, k = 10, fill = "extend")
Lam_hi <- rollmean(precs_hi, k = 10, fill = "extend")
Lam_pt <- rollmean(precs_pt, k = 10, fill = "extend")
Lam_z <- Lam_pt + Lam_lo + Lam_hi

# Get responses
y_pt <- pk2[, 2]
y_lo <- as.matrix(pk2[, 3:18]) # each column is one run
y_hi <- pk2[, 19]

# Transform to emulation space
y_pt <- scrP(y_pt, k)
y_lo <- scrP(y_lo, k)
y_hi <- scrP(y_hi, k)

# Get averages
y_lra <- rowMeans(y_lo) # low resolution average
y_avg <- (1 / Lam_z) * (Lam_pt * y_pt + Lam_lo * y_lra + Lam_hi * y_hi)

# Load in fitted values
fit <- read.csv(paste0("../fitting/results/dgp_", model_name, ".csv"))

###################
# Plot data types #
###################

# color-blind friendly color palette
cbcols <- palette.colors(palette = "Okabe-Ito")
#1 black, 2 orange, 3 skyblue , 4 bluishgreen, 5 yellow 
#6 blue, 7 vermillion,  8 reddishpurple, 9 gray 

# plot_data.png
if(JPG) jpeg("../paper/plot_data.jpeg", width = 12, height = 4, units = "in", res = 300)
par(mfrow=c(1,3), mar=c( 5.6, 4.6, 4.6, 2.1))
matplot(log10(k), y_lo, col = cbcols[9], 
        type = "l", lty = 1, ylim = range(y_avg), lwd=0.5,
        xlab=expression(log[10](k)), cex.axis=1.75, cex.lab=1.75,
        ylab=expression(paste("\U1D4AB","  (k)")))
lines(log10(k), y_pt, col = cbcols[4])
lines(log10(k), y_hi, col = cbcols[8], lwd=2)
lines(log10(k)[index_list$pert.ix], rep(-1, length(index_list$pert.ix)), 
      col = cbcols[4], lty = 3, lwd=3)
lines(log10(k)[index_list$lowres.ix], rep(-1.1, length(index_list$lowres.ix)), 
      col = cbcols[9], lty = 3, lwd=3)
lines(log10(k)[index_list$highres.ix], rep(-1.2, length(index_list$highres.ix)), 
      col = cbcols[8], lty = 3, lwd=3)
lines(log10(k), y_avg, col = cbcols[2], lty = 2, lwd = 3)

# Same as first plot, but only for the low-resolution indices
matplot(log10(k)[index_list$lowres.ix], y_lo[index_list$lowres.ix,], col = cbcols[9], 
        type = "l", lty = 1, ylim = range(y_avg[index_list$lowres.ix])+c(-.05,.05), lwd=0.5,
        xlab=expression(log[10](k)), cex.axis=1.75, cex.lab=1.75,
        ylab=expression(paste("\U1D4AB","  (k)")))
lines(log10(k), y_pt, col = cbcols[4])
lines(log10(k)[index_list$lowres.ix], y_hi[index_list$lowres.ix], col = cbcols[8], lwd=2)
lines(log10(k)[index_list$lowres.ix], y_avg[index_list$lowres.ix], 
      col = cbcols[2], lty = 2, lwd = 3)
legend(x = "topleft", legend = c("pert","low res","hi res", "wt avg"),
       col = cbcols[c(4,9,8,2)], lty = c(1,1,1,2), 
       lwd = c(1,0.5,2,3), cex=1.25,
       inset = 0.01, bty="n"  # nudges legend inward
)

loess_fit <- loess(y_avg ~ x, span = 0.15)

# same as first plot, but LOESS-smoothed weighted average subtracted
matplot(log10(k), y_lo - loess_fit$fitted, col = cbcols[9], xlim = c(log10(0.04),log10(0.25)),
        type = "l", lty = 1, ylim = c(-.02, .02), lwd=0.5,
        xlab=expression(log[10](k)), cex.axis=1.75, cex.lab=1.75,
        ylab=expression(paste("\U1D4AB","  (k)", ", centered")))
lines(log10(k), y_pt-loess_fit$fitted, col = cbcols[4])
lines(log10(k), y_hi-loess_fit$fitted, col = cbcols[8], lwd=2)
lines(log10(k), y_avg-loess_fit$fitted, 
      col = cbcols[2], lty = 2, lwd = 3)
if(JPG) dev.off()

######################
# Plot fitted values #
######################

par(mfrow=c(1,1))
# Plot model fit alongside data
matplot(log10(k), y_lo, type="l", lty=1,
        col=cbcols[9], xlab=expression(log[10](k)), 
        ylab=expression(paste("\U1D4AB","(k)")),
        ylim=range(fit$ub, fit$lb))
lines(log10(k),y_pt, col=cbcols[4], lwd=2, lty=3)
lines(log10(k),y_hi, col=cbcols[8], lwd=2, lty=3)
lines(log10(k), y_avg, col=cbcols[2], lwd=2, lty=2)
lines(log10(k),fit$ub, col=cbcols[1],lty=2,lwd=2)
lines(log10(k),fit$lb, col=cbcols[1],lty=2,lwd=2)
legend(x = "bottomright", legend = c("pert","low res","hi res", "wt avg", "UQ"),
       col = cbcols[c(4,9,8,2,1)], lty = c(3,1,3,2,2,2), 
       lwd = c(2,1,2,2,2))

# plot_fit.png
if(JPG) jpeg("../paper/plot_fit.jpeg", width = 12, height = 4, units = "in", res = 300)

# Set layout: 2 columns, relative widths 2:1
layout(matrix(1:2, nrow = 1), widths = c(1, 1))

# Left plot (2/3 width)
par(mar = c(4, 4, 2, 1))  # set margins: bottom, left, top, right

# Plot model fit alongside data (posterior mean removed)
matplot(log10(k), y_lo - fit$m, type="l", lty=3,
        col=cbcols[9], ylim = c(-.02,.02), xlab=expression(log[10](k)), 
        ylab=expression(paste("\U1D4AB","  (k), mean removed")))
abline(h=0, col=cbcols[1],lty=1,lwd=1)
lines(log10(k), y_pt - fit$m, col=cbcols[4], lwd=3, lty=4)
lines(log10(k), y_hi - fit$m, col=cbcols[8], lwd=3, lty=2)
lines(log10(k), y_avg - fit$m , col=cbcols[2], lwd=2, lty=1)
lines(log10(k), fit$ub - fit$m, col=cbcols[1],lty=2,lwd=2)
lines(log10(k), fit$lb - fit$m, col=cbcols[1],lty=2,lwd=2)
legend(x = "bottomright", legend = c("pert","low res","hi res", "wt avg", "UQ"),
       col = cbcols[c(4,9,8,2,1)], lty = c(4,3,2,1,2), 
       lwd = c(2,1,2,2,2), cex=0.9, bty="n")


# Right plot (1/3 width) - keep square aspect

par(mar = c(4, 4, 2, 1))
# load results from burned in fit for M001
load("../fitting/results/fit_M001.rda")

indx <- floor(seq(from = 1, to = fit$nmcmc, length = 100))
if (indx[1] == 0) indx[1] <- 1
o <- order(fit$x)
w <- fit$w[o, indx]
for (i in 1:length(indx)) w[, i] <- w[, i] - mean(w[, i])
matplot(log10(k), w, type = "l", xlab = expression(paste(log[10](k))), ylab = "W", col = "grey", 
        ylim = c(min(w), max(w)))
if(JPG) dev.off()


