
# Plot One Layer --------------------------------------------------------------
#' @rdname plot
#' @export

plot.gphm <- function(x) {
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  Dx <- ncol(x$x)
  
  # Trace plots
  par(mfrow = c(1, 2), mar = c(5, 4, 2, 2))
  plot(x$ll, type = 'l', ylab = 'logl', xlab = 'Iteration',
       main = 'Trace Plot of logl')
  plot(x$theta, type = 'l', ylab = 'theta_y', xlab = 'Iteration',
      main = 'Trace Plot of theta')
      
}

# Plot Two Layer --------------------------------------------------------------
#' @rdname plot
#' @export

plot.dgp2hm <- function(x) {
  
  # save and restore par settings 
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  Dx <- ncol(x$x)
  D <- ncol(x$w[[1]])
  
  # Trace plots
  par(mfrow = c(1, D + 2), mar = c(5, 4, 2, 2))
  plot(x$ll, type = 'l', ylab = 'outer logl', xlab = 'Iteration',
       main = 'Trace Plot of outer logl')
  plot(x$theta_y, type = 'l', ylab = 'theta_y', xlab = 'Iteration',
       main = 'Trace Plot of theta_y')
  for (i in 1:D)
    plot(x$theta_w[, i], type = 'l', ylab = 'theta_w', xlab = 'Iteration',
         main = paste0('Trace Plot of theta_w [', i, ']'))
  
  # Hidden layer plot
  indx <- floor(seq(from = 1, to = x$nmcmc, length = 100))
  if (indx[1] == 0) indx[1] <- 1
  col <- heat.colors(100 + 10) # add ten to avoid using colors that are too light
  par(mfrow = c(1, 1), mar = c(4, 4, 2, 2))
  o <- order(x$x)
  plot(x$x[o], x$w[[indx[1]]][o] - mean(x$w[[indx[1]]]), type = 'l', xlab = 'X', 
       ylab = 'W', col = col[1], main = paste0('MCMC samples of X to W'), 
       ylim = c(min(unlist(x$w[indx])), max(unlist(x$w[indx]))))
  for (j in 2:length(indx)) {
    lines(x$x[o], x$w[[indx[j]]][o] - mean(x$w[[indx[j]]]), col = col[j])
  }

}
