
# Function Contents -----------------------------------------------------------
# External:
#   trim.gphm
#   trim.dgp2hm 

#' @rdname trim
#' @export
trim <- function(object, burn, thin)
  UseMethod("trim", object)

# Trim One Layer --------------------------------------------------------------
#' @rdname trim
#' @export

trim.gphm <- function(object, burn, thin = 1) {
  
  tic <- proc.time()[[3]]
  
  if (burn >= object$nmcmc) stop("burn must be less than nmcmc")
  
  nmcmc <- object$nmcmc
  indx <- (burn + 1):nmcmc
  indx <- indx[which(indx %% thin == 0)]
  
  object$nmcmc <- length(indx)
  object$theta <- object$theta[indx, drop = FALSE]
  object$tau2 <- object$tau2[indx, drop = FALSE]
  object$ll <- object$ll[indx, drop = FALSE]
  
  toc <- proc.time()[[3]]
  object$time <- object$time + unname(toc - tic)
  
  return(object)
}

# Trim Two Layer --------------------------------------------------------------
#' @rdname trim
#' @export

trim.dgp2hm <- function(object, burn, thin = 1) {
  
  tic <- proc.time()[[3]]
  
  if (burn >= object$nmcmc) stop("burn must be less than nmcmc")
  
  nmcmc <- object$nmcmc
  indx <- (burn + 1):nmcmc
  indx <- indx[which(indx %% thin == 0)]
  
  object$nmcmc <- length(indx)
  object$theta_y <- object$theta_y[indx, drop = FALSE]
  object$theta_w <- object$theta_w[indx, , drop = FALSE]
  object$w <- as.list(object$w[indx])
  object$tau2 <- object$tau2[indx, drop = FALSE]
  object$ll <- object$ll[indx, drop = FALSE]
  
  toc <- proc.time()[[3]]
  object$time <- object$time + unname(toc - tic)
  
  return(object)
}