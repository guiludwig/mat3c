#' Predictive Analysis for a fitmat3 object
#'
#' Evaluates the predictive posterior probability of T > T_obs for various statistics.
#'
#' @param model An fitmat3 object (resulting from a call to fitmat3)
#' @param fname If the MCMC samples were saved into text files, this can be used to generate
#'              predictive probabilities from them. Defaults to NULL.
#' @param whichSamples Which (and how many) samples to collect in order to generate the
#'              predictive probabilities. Must be smaller than the number of MCMC samples.
#' @param additionalArgs additional arguments to be passed to rmat3() (if read from text files,
#'              it is necessary to provide a list with window, R and R_c).
#'
#' @export
#' @return \code{predictiveStats} returns a list containing at least the following components
#'   \item{pL}{predictive probability for the expected number of centers}
#'   \item{pRc}{predictive probability for the number of fingers within R_centers distance of other centers}
#'   \item{pnx}{predictive probability for the average number of fingers}
#'   \item{psx}{predictive probability for the standard deviation of the number of fingers}
#'   \item{chx}{predictive probability for the mean area of the convex hull of fingers (reactive territory)}
#'   \item{spx}{predictive probability for the spanned angle of fingers}
#'   \item{mdx}{predictive probability for the average of minimum distance between fingers}
#'   \item{adx}{predictive probability for the average of average distance between fingers}
#'   \item{cex}{predictive probability for the average distance between centers}
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
  mchull <- mean(sapply(fmat.cond, function(x) areahull(cbind(x$x, x$y))))
  spangle <- mean(sapply(fmat.cond, spanAngle))
  minD <- mean(sapply(fmat.cond, minDist))
  meanD <- mean(sapply(fmat.cond, meanDist))
  meanCD <- mean(dist(t(sapply(fmat.cond, function(x) x$centers, simplify = "matrix"))))
  if(length(meanCD) == 0) {
    meanCD <- 0
  }

  M <- length(whichSamples)
  suff0 <- matrix(0, ncol = 2, nrow = M)
  n.fingers0 <- numeric(M)
  s.fingers0 <- numeric(M)
  mchull0 <- numeric(M)
  spangle0 <- numeric(M)
  minD0 <- numeric(M)
  meanD0 <- numeric(M)
  meanCD0 <- numeric(M)

  for(i in seq_along(whichSamples)){
    cat(sprintf("%3.1f\n", 100*i/M))
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
    mchull0[i] <- mean(sapply(x, function(x) areahull(x$fingers)))
    spangle0[i] <- mean(sapply(x, spanAngle))
    minD0[i] <- mean(sapply(x, minDist))
    meanD0[i] <- mean(sapply(x, meanDist))
    tempCD <- mean(dist(t(sapply(x, function(x) x$centers, simplify = "matrix"))))
    if(length(tempCD) == 0) {
      tempCD <- 0
    }
    meanCD0[i] <- tempCD
  }
  suff0[,2] <- -1*suff0[,2]
  psuff1 <- sum(suff[1] < suff0[,1])/M
  psuff2 <- sum(suff[2] < suff0[,2])/M
  pn <- sum(n.fingers < n.fingers0)/M
  ps <- sum(s.fingers < s.fingers0)/M
  ch <- sum(mchull < mchull0)/M
  sp <- sum(spangle < spangle0)/M
  md <- sum(minD < minD0)/M
  ad <- sum(meanD < meanD0)/M
  cd <- sum(meanCD < meanCD0)/M

  return(list(pL = psuff1,
              pRc = psuff2,
              pnx = pn,
              psx = ps,
              chx = ch,
              spx = sp,
              mdx = md,
              adx = ad,
              cex = cd))
}
