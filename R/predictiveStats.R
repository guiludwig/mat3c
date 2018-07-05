#' @export
predictiveStats <- function(model, whichSamples = seq(1, nrow(model$parameters), 10)){
  fmat.cond <- model$augmentedData
  suff <- sufficientStat(model$augmentedData, model$R_centers)
  suff[2] <- -1*suff[2]
  n.fingers <- mean(sapply(model$augmentedData, function(x) length(x$x)))
  s.fingers <- sd(sapply(model$augmentedData, function(x) length(x$x)))
  L <- length(whichSamples)
  suff0 <- matrix(0, ncol = 2, nrow = L)
  n.fingers0 <- numeric(L)
  s.fingers0 <- numeric(L)
  for(i in seq_along(whichSamples)){
    cat(sprintf("%3.1f\n", 100*i/L))
    x <- rmat3(beta  = model$parameters[whichSamples[i], "beta"],
               phi   = model$parameters[whichSamples[i], "phi"],
               gamma = model$parameters[whichSamples[i], "gamma"],
               kappa = model$parameters[whichSamples[i], "kappa"],
               sigma = model$parameters[whichSamples[i], "sigma"],
               win = model$win,
               R_centers = model$R_centers,
               R_clusters = model$R_clusters)
    suff0[i, ] <- sufficientStat(x, model$R_centers)
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
