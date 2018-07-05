#' Plot likelihood function
#'
#' Plots the likelihood function for a mat3 object
#'
#' @param final ...
#' @param psi ...
#' @param win ...
#' @param R_clusters ...
#' @param R_centers ...
#' @param J ...
#' @param N ...
#' @param logPriors ...
#' @param initialValues ...
#' @param GeyerThompson ...
#' @param fname ...
#' @param logv ...
#'
#' @export
#' @examples
#' \notrun{x <- rmat3(70, 2, 5, 0.05, 3)
#' gr <- plotLikelihood(x, GeyerThompson = "importance")
#' # No importance sampling:
#' gr <- plotLikelihood(x, psi = c(70, 1, 5, 0.05, 3), GeyerThompson = "simple")
#' # Reading from a file
#' gr <- plotLikelihood(NA, fname = c("centers-param=3-b=19.txt", "fingers-param=3-b=19.txt"))}
plotLikelihood <- function(final,
                           psi = c(beta = 70,
                                   phi = 2,
                                   gamma = 5,
                                   sigma = 0.05,
                                   kappa = 3),
                           win = owin(c(0,1), c(0,1)),
                           R_clusters = 0.02, R_centers = 0.05,
                           J = 100, N = 10000,
                           logPriors = preparePriors(),
                           initialValues = pickInitialValues(),
                           GeyerThompson = c("large", "importance", "simple"),
                           fname = NULL,
                           logv = TRUE){

  GeyerThompson <- match.arg(GeyerThompson)

  if(!is.na(final[1])){
    fmat <- lapply(final, function(x) list(centers = x$centers, x = x$fingers[,1], y = x$fingers[,2]))
  } else if(!is.null(fname)) {
    centers <- read.table(fname[1])
    fingers <- read.table(fname[2])
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
    stop("Provide data (from rmat3 object) or a link to centers/fingers.")
  }

  # Preparation

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

  theta.fixo <- log(c(beta.p, phi.p))

  suff.mcmc <- switch(GeyerThompson,
                      large = if(as.numeric(psi[2]) == 1){
                        sampleZ1s(psi[1], psi[2], psi[3],
                                  psi[4], psi[5],
                                  win, J*100, R_centers, R_clusters,
                                  debug = TRUE)
                      } else {
                        sampleZ1(psi[1], psi[2], psi[3],
                                 psi[4], psi[5],
                                 win, J*100, R_centers, R_clusters,
                                 debug = TRUE)
                      },
                      importance = if(as.numeric(psi[2]) == 1){
                        sampleZ1s(psi[1], psi[2], psi[3],
                                  psi[4], psi[5],
                                  win, J, R_centers, R_clusters,
                                  debug = TRUE)
                      } else {
                        sampleZ1(psi[1], psi[2], psi[3],
                                 psi[4], psi[5],
                                 win, J, R_centers, R_clusters,
                                 debug = TRUE)
                      },
                      simple = sampleZ1(psi[1], psi[2], psi[3],
                                        psi[4], psi[5],
                                        win, J, R_centers, R_clusters,
                                        debug = TRUE)
  )

  fmat.cond <- sampleThin(fmat, psi[3], psi[5],
                          psi[4], win, R_clusters, up.0)

  # png("estimationZ1.png")
  layout(matrix(1:2, ncol=2))
  barplot(table(suff.mcmc[,1]), main = "suff.stat1")
  abline(v = which(suff.beta.phi[1] == as.numeric(names(table(suff.mcmc[,1])))), col="Red")
  barplot(table(suff.mcmc[,2]), main = "suff.stat2")
  abline(v = which(suff.beta.phi[2] == as.numeric(names(table(suff.mcmc[,2])))), col="Red")
  # dev.off()

  #!# DEBUG FROM HERE #!#

  f_beta <- function(x) posterioriBeta(fmat.cond, beta = x,
                                       phi = phi.p, gamma = gamma.p,
                                       sigma = sigma.p, kappa = kappa.p,
                                       BETA.P = beta.p,
                                       SUFF.BETA.PHI = suff.beta.phi,
                                       THETA.FIXO = theta.fixo,
                                       SUFF.MCMC = suff.mcmc,
                                       SFALL = suff.all,
                                       logpriorBeta = logPriors$priorBeta)
  f_phi <- function(x) posterioriPhi(fmat.cond, beta = beta.p,
                                     phi = x, gamma = gamma.p,
                                     sigma = sigma.p, kapp = kappa.p,
                                     SUFF.BETA.PHI = suff.beta.phi,
                                     THETA.FIXO = theta.fixo,
                                     SUFF.MCMC = suff.mcmc,
                                     SFALL = suff.all,
                                     logpriorPhi = logPriors$priorPhi)
  f_gamma <- function(x) posterioriGamma(fmat.cond, beta = beta.p,
                                         phi = phi.p, gamma = x,
                                         sigma = sigma.p, kappa = kappa.p,
                                         GAMMA.P = gamma.p,
                                         SFALL = suff.all,
                                         logpriorGamma = logPriors$priorGamma)
  f_sigma <- function(x) posterioriSigma(fmat.cond, beta = beta.p,
                                         phi = phi.p, gamma = gamma.p,
                                         sigma = x, kappa = kappa.p,
                                         SIGMA.P = sigma.p,
                                         A1 = a1, A2 = a2, B1 = b1, B2 = b2,
                                         NN = N, Rc = R_centers, R = R_clusters,
                                         logpriorSigma = logPriors$priorSigma)
  f_kappa <- function(x) posterioriKappa(fmat.cond, beta = beta.p,
                                         phi = phi.p, gamma = gamma.p,
                                         sigma = sigma.p, kappa = x,
                                         KAPPA.P = kappa.p,
                                         A1 = a1, A2 = a2, B1 = b1, B2 = b2,
                                         NN = N, Rc = R_centers, R = R_clusters,
                                         logpriorKappa = logPriors$priorKappa)

  I <- 200
  beta0  <- seq(0.1*psi[1], 10*psi[1], length.out = I) # 70
  phi0   <- seq(0.1*psi[2], 10*psi[2], length.out = I) # 2
  gamma0 <- seq(0.1*psi[3], 10*psi[3], length.out = I) # 50
  sigma0 <- seq(0.1*psi[4], 10*psi[4], length.out = I) # 0.025
  kappa0 <- seq(0.1*psi[5], 10*psi[5], length.out = I) # 3

  betaPos <- betaPrior <- betaLik <- numeric(I)
  phiPos <- phiPrior <- phiLik <- numeric(I)
  gammaPos <- gammaPrior <- gammaLik <- numeric(I)
  sigmaPos <- sigmaPrior <- sigmaLik <- numeric(I)
  kappaPos <- kappaPrior <- kappaLik <- numeric(I)

  suff.all <- sufficientStatAll(fmat.cond, sigma.p, kappa.p,
                                A1 = a1, A2 = a2, B1 = b1, B2 = b2,
                                NN = N, Rc = R_centers, R = R_clusters)

  for(i in seq_len(I)){
    print(i)
    # BETA
    betaPos[i] <- f_beta(beta0[i])
    betaPrior[i] <- preparePriors()$priorBeta(beta0[i], beta.p)
    betaLik[i] <- betaPos[i] - betaPrior[i]
    # PHI
    phiPos[i] <- f_phi(phi0[i])
    phiPrior[i] <- preparePriors()$priorPhi(phi0[i])
    phiLik[i] <- phiPos[i] - phiPrior[i]
    # GAMMA
    gammaPos[i] <- f_gamma(gamma0[i])
    gammaPrior[i] <- preparePriors()$priorGamma(gamma0[i], gamma.p)
    gammaLik[i] <- gammaPos[i] - gammaPrior[i]
    # SIGMA
    sigmaPos[i] <- f_sigma(sigma0[i])
    sigmaPrior[i] <- preparePriors()$priorSigma(sigma0[i], sigma.p)
    sigmaLik[i] <- sigmaPos[i] - sigmaPrior[i]
    # KAPPA
    kappaPos[i] <- f_kappa(kappa0[i])
    kappaPrior[i] <- preparePriors()$priorKappa(kappa0[i], kappa.p)
    kappaLik[i] <- kappaPos[i] - kappaPrior[i]
  }

  # Actual Plot

  if(logv) {
    f <- function(x) x
  } else {
    f <- function(x) exp(x)
  }

  layout(matrix(c(1:6), ncol=2))
  # BETA
  plot(beta0, f(betaPos), type = "l", ylim = range(f(c(betaPos,betaPrior,betaLik))))
  lines(beta0, f(betaLik), col="Blue")
  abline(v = psi[1], col="Blue", lty = 2)
  lines(beta0, f(betaPrior), col="Red")
  abline(v = beta.p, col="Red", lty = 2)
  # PHI
  plot(phi0, f(phiPos), type = "l", ylim = range(f(c(phiPos,phiPrior,phiLik))))
  lines(phi0, f(phiLik), col="Blue")
  abline(v = psi[2], col="Blue", lty = 2)
  lines(phi0, f(phiPrior), col="Red")
  abline(v = phi.p, col="Red", lty = 2)
  # GAMMA
  plot(gamma0, f(gammaPos), type = "l", ylim = range(f(c(gammaPos,gammaPrior,gammaLik))))
  lines(gamma0, f(gammaLik), col="Blue")
  abline(v = psi[3], col="Blue", lty = 2)
  lines(gamma0, f(gammaPrior), col="Red")
  abline(v = gamma.p, col="Red", lty = 2)
  # SIGMA
  plot(sigma0, f(sigmaPos), type = "l", ylim = range(f(c(sigmaPos,sigmaPrior,sigmaLik))))
  lines(sigma0, f(sigmaLik), col="Blue")
  abline(v = psi[4], col="Blue", lty = 2)
  lines(sigma0, f(sigmaPrior), col="Red")
  abline(v = sigma.p, col="Red", lty = 2)
  # KAPPA
  plot(kappa0, f(kappaPos), type = "l", ylim = range(f(c(kappaPos,kappaPrior,kappaLik))))
  lines(kappa0, f(kappaLik), col="Blue")
  abline(v = psi[5], col="Blue", lty = 2)
  lines(kappa0, f(kappaPrior), col="Red")
  abline(v = kappa.p, col="Red", lty = 2)
  # LEGEND
  plot(0, type="n", axes=FALSE, xlab="", ylab="")
  legend("center", lty=c(1,1,1,2,2), col=c("Black", "Blue", "Red","Blue", "Red"),
         legend = c("Posterior", "Likelihood", "Prior", "Truth", "Starting value"))

  invisible(list(beta = list(beta0, betaPos, betaLik),
                 phi = list(phi0, phiPos, phiLik),
                 gamma = list(gamma0, gammaPos, gammaLik),
                 sigma = list(sigma0, sigmaPos, sigmaLik),
                 kappa = list(kappa0, kappaPos, kappaLik)))

}
