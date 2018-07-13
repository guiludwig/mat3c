#' @export
predictiveStats <- function(model, fname = NULL, whichSamples, additionalArgs = list(win = owin(),
                                                                                      R = 0.005,
                                                                                      R_c = 0.02), ...){
  if(inherits(model, "fitmat3")){
    fmat.cond <- model$augmentedData
    parameters <- model$parameters
    R_c <- model$R_centers
    R <- model$R_clusters
    win <- model$win
  } else {
    if(is.null(fname)) stop("If model isn't of class fitmat3, must provide formatted input")
    centers <- read.table(fname[1], ...)
    fingers <- read.table(fname[2], ...)
    parameters <- read.table(fname[3], ...)
    colnames(parameters) <- c("beta", "phi", "gamma", "sigma", "kappa")
    fmat <- vector("list", length(centers$Tree))
    final <- vector("list", length(centers$Tree))
    for (i in 1:length(centers$Tree)) {
      temp <- matrix(cbind(fingers$X[fingers$Tree == i],
                           fingers$Y[fingers$Tree == i]),
                     ncol=2)
      final[[i]] <- list(centers = c(centers$X[i],centers$Y[i]),
                         fingers = temp)
      fmat[[i]] <- list(centers = c(centers$X[i],centers$Y[i]),
                        x = temp[,1], y = temp[,2])
    }
    fmat.cond <- fmat
    R_c <- additionalArgs$R_c
    R <- additionalArgs$R
    win <- additionalArgs$win
  }

  suff <- sufficientStat(fmat.cond, R_c)
  suff[2] <- -1*suff[2]
  n.fingers <- mean(sapply(fmat.cond, function(x) length(x$x)))
  s.fingers <- sd(sapply(fmat.cond, function(x) length(x$x)))
  L <- length(whichSamples)
  suff0 <- matrix(0, ncol = 2, nrow = L)
  n.fingers0 <- numeric(L)
  s.fingers0 <- numeric(L)
  for(i in seq_along(whichSamples)){
    cat(sprintf("%3.1f\n", 100*i/L))
    x <- rmat3(beta  = parameters[whichSamples[i], "beta"],
               phi   = parameters[whichSamples[i], "phi"],
               gamma = parameters[whichSamples[i], "gamma"],
               kappa = parameters[whichSamples[i], "kappa"],
               sigma = parameters[whichSamples[i], "sigma"],
               win = win,
               R_centers = R_c,
               R_clusters = R)
    suff0[i, ] <- sufficientStat(x, R_c)
    n.fingers0[i] <- mean(sapply(x, function(x) nrow(x$fingers)))
    s.fingers0[i] <- sd(sapply(x, function(x) nrow(x$fingers)))
  }
  suff0[,2] <- -1*suff0[,2]
  psuff1 <- sum(suff[1] < suff0[,1])/L
  psuff2 <- sum(suff[2] < suff0[,2])/L
  pn <- sum(n.fingers < n.fingers0)/L
  ps <- sum(s.fingers < s.fingers0)/L
  return(list(pL = psuff1,
              pRc = psuff2,
              pnx = pn,
              psx = ps))
}
