#' Fit a Matern-III marked point process model via MCMC sampling
#'
#' Markov chain Monte Carlo sampler from Matern-III marked spatial point process
#'
#' @param mat3 A mat3 object, see \code{\link{rmat3}}. If reading from a data file,
#'              set to NA instead.
#' @param fname Will only be read if mat3 is set to NA. File names with the centers
#'              and fingers dataset, to be read with \code{\link{read.table}}. The default
#'              action assumes the files "centers.txt" and "fingers.txt", each formatted
#'              with a header and three columns corresponding to a center ID (integer),
#'              and the spatial location of the center/finger end. The fingers file should
#'              have at least one finger for each center ID in the centers file.
#' @param win Square window in which the process is to be sampled. See \code{\link{owin}}.
#'              Defaults to the unit square.
#' @param R_clusters Nuisance parameter for the fingers' inhibition. No finger endings
#'              will be closer than R_clusters in Euclidean distance during the process
#'              sampling step (hard core inhibiition).
#' @param R_centers Nuisance parameter for the centers' inhibition. At the birth-and-death
#'              process, process will take into account how many fingers (attached to other
#'              active centers) are currently at R_centers's distance from the candidate
#'              center being born.
#' @param N Number of samples for integrating the area of the shadow. Defaults to 10000.
#' @param seed Fixes the RNG seed for replication. Defaults to NULL, which does not
#'              fix the seed.
#' @param resultsName ...
#' @param initialValues ...
#' @param ... further arguments passed to \code{\link{read.table}}.
#'
#' @export
#' @useDynLib mat3c
#' @importFrom Rcpp evalCpp
#' @return \code{fitmat3} returns a list containing at least the following components
#'   \item{parameters}{An L by 5 matrix with the samples from beta, phi, gamma, sigma, kappa.}
#'
#' @examples
#' set.seed(1234)
#' x <- rmat3(70, 2, 5, 0.05, 3)
#' plot(x)
#' # Changing default sampling sizes to make it run fast
#' model <- pseudofitmat3(x, N = 1000, seed = 1234)
#'
#' @author Guilherme Ludwig and Nancy Garcia
#'
#' @references
#'
#'   Garcia, N., Guttorp, P. and Ludwig, G. (2018) TBD
#'
#' @seealso \code{\link{rmat3}}
#' @keywords Spatial Statistics
#' @keywords Functional Data Analysis
pseudofitmat3 <- function(mat3, fname = c("centers.txt", "fingers.txt"),
                    win = owin(c(0,1), c(0,1)),
                    R_clusters = 0.005, R_centers = 0.02,
                    N = 10000, seed = NULL,
                    resultsName = NULL,
                    initialValues = pickInitialValues(),
                    ...){

  set.seed(seed)

  if(is.na(mat3[1])){
    centers <- read.table(fname[1], ...)
    fingers <- read.table(fname[2], ...)
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
  } else {
    final <- mat3
    fmat <- lapply(mat3, function(x) list(centers = x$centers, x = x$fingers[,1], y = x$fingers[,2]))
  }
  a1 <- win$x[1]
  a2 <- win$x[2]
  b1 <- win$y[1]
  b2 <- win$y[2]
  limx <- c(a1, a2)
  limy <- c(b1, b2)

  suff.beta.phi <- sufficientStat(final, Rc = R_centers)

  Ln <- length(fmat)
  n.j <- sapply(fmat, function(x) length(x$x))

  up.0 <- rep(0, Ln)
  angles <- NULL
  for (ij in 1:length(fmat)){
    temp <- atan2(fmat[[ij]]$y-fmat[[ij]]$centers[2],
                  fmat[[ij]]$x-fmat[[ij]]$centers[1])
    temp[temp < 0] <- 2*pi + temp[temp < 0]
    up.0[ij] <- mean(circular(temp, modulo = "2pi"))
    angles <- c(angles, circular(temp - up.0[ij])) # Non-circular, mean 0
  }

  gamma.p <- initialValues$iniGamma(fmat)
  beta.p  <- initialValues$iniBeta(fmat, win)
  phi.p   <- initialValues$iniPhi()
  kappa.p <- initialValues$iniKappa(angles)
  sigma.p <- initialValues$iniSigma(fmat)

  fmat.ppp <- as.ppp(t(sapply(fmat, function(x) x$centers)), win)
  Q <- quadscheme(fmat.ppp)

  # I need:
  # evaluate Papangelou intensity at dummy points
  # pagangelouIntensity(beta, phi, Q, fmat, R_centers, Boot)




  return(list(parameters = parameters,
              window = win,
              R_centers = R_centers,
              R_clusters = R_clusters))

}
