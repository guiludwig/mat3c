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
#' @param L Number of samples for the MCMC sampler. Defaults to 10000.
#' @param N Number of samples for integrating the area of the shadow. Defaults to 10000.
#' @param J Number of updates for the Geyer-Thompson step in estimating the sampling
#'              distribution for beta and phi (if GeyerThompson is set to "large", 100 times J).
#'              Defaults to 100 for "simple" and "importance" options in GeyerThompson,
#'              and 10000 for "large".
#' @param seed Fixes the RNG seed for replication. Defaults to NULL, which does not
#'              fix the seed.
#' @param resultsName Writes the results to a file; defaults to NULL, which does not
#'              save the results.
#' @param logPriors A list of 5 log-priors for the coefficients, each being a function of the
#'              parameter and initial value only. See \code{\link{preparePriors}}, which
#'              produces the default values of logPriors.
#' @param initialValues Some holistic initial values for empirical priors on the parameters.
#'              See \code{\link{pickInitialValues}}, which produces the default values of
#'              initialValues.
#' @param GeyerThompson Provides different methods for updating the Geyer and Thompson (1999)
#'              step in the computation of the likelihood. The "large" option produces J*100
#'              samples from which the Geyer-Thompson approximation ot the likelihood is
#'              built; "importance" uses J samples, but they are updated at every 100
#'              iterations. "simple" works like "importance", but sets phi = 1. "Resample"
#'              uses a resample method on existing finger. "MollerAlt" uses an alternative
#'              approach of auxiliary variabiles to avoid computing the normalizing constant,
#'              described in Moller et al (2006).
#' @param tildeTheta Only needed for the "MollerAlt" option in "GeyerThompson". Corresponds
#'              to a vector of length two, with estimates to parameters beta and phi.
#'              Defaults to average intensity and 1.
#' @param candidateVar Variability of MCMC candidate selection step, if tuning is necessary.
#' @param verbose Prints acceptance rates for candidates. Defaults to TRUE.
#' @param ... further arguments passed to \code{\link{read.table}}.
#'
#' @export
#' @useDynLib mat3c
#' @import circular spatstat
#' @importFrom Rcpp evalCpp
#' @return \code{fitmat3} returns a list containing at least the following components
#'   \item{parameters}{An L by 5 matrix with the samples from beta, phi, gamma, sigma, kappa.}
#'
#' @examples
#' set.seed(1234)
#' x <- rmat3(70, 2, 5, 0.05, 3)
#' plot(x)
#' # Changing default sampling sizes to make it run fast
#' model <- fitmat3(x, L = 100, N = 1000, J = 2, seed = 1234)
#' burnin <- 20
#' layout(matrix(c(1:5,0), ncol=2))
#' plot(model$parameters[-c(1:burnin),1], type="l", ylab="beta")
#' abline(h = median(model$parameters[-c(1:burnin),1]))
#' abline(h = 2.5*25, col = "Red")
#' plot(model$parameters[-c(1:burnin),2], type="l", ylab="phi")
#' abline(h = median(model$parameters[-c(1:burnin),2]))
#' abline(h = 2, col = "Red")
#' plot(model$parameters[-c(1:burnin),3], type="l", ylab="gamma")
#' abline(h = median(model$parameters[-c(1:burnin),3]))
#' abline(h = 20, col = "Red")
#' plot(model$parameters[-c(1:burnin),4], type="l", ylab="sigma")
#' abline(h = median(model$parameters[-c(1:burnin),4]))
#' abline(h = 0.05, col = "Red")
#' plot(model$parameters[-c(1:burnin),5], type="l", ylab="kappa")
#' abline(h = median(model$parameters[-c(1:burnin),5]))
#' abline(h = 3, col = "Red")
#' # Saving results to file:
#' \dontrun{model <- fitmat3(NA, fname = c("centers.txt", "fingers.txt"),
#'                           seed = 1234, resultsName = "results.txt")}
#'
#' @author Guilherme Ludwig and Nancy Garcia
#'
#' @references
#'
#'   Garcia, N., Guttorp, P. and Ludwig, G. (2020) "Interacting cluster point process model
#'   for epidermal nerve fibers", to appear.
#'
#'   Geyer, C. J., and Thompson, E. A. (1992) "Constrained Monte Carlo maximum likelihood
#'   for dependent data." Journal of the Royal Statistical Society, Series B, Vol. 54,
#'   pp. 657-699.
#'
#'   Moller, J., Pettitt, A. N., Reeves, R. and Berthelsen, K. K. (2006) "An efficient
#'   Markov chain Monte Carlo method for distributions with intractable normalising
#'   constants." Biometrika, Vol. 93, No. 2, pp. 451-458.
#'
#' @seealso \code{\link{rmat3}}
#' @keywords Spatial statistics
#' @keywords Spatial point processes
#' @keywords Bayesian methods
fitmat3 <- function(mat3, fname = c("centers.txt", "fingers.txt"),
                    win = owin(c(0,1), c(0,1)),
                    R_clusters = 0.005, R_centers = 0.02,
                    L = 10000, N = 10000, J = 100, seed = NULL,
                    resultsName = NULL,
                    logPriors = preparePriors(),
                    initialValues = pickInitialValues(),
                    GeyerThompson = c("large", "resample",
                                      "importance", "simple", "MollerAlt"),
                    tildeTheta = c(length(mat3)/area(win), 1),
                    candidateVar = c(0.1, 0.2, 0.1, 0.1, 0.1),
                    verbose = TRUE,
                    ...){

  set.seed(seed)

  priornames <- c("priorBeta", "priorPhi", "priorGamma", "priorSigma", "priorKappa")
  stopifnot(all.equal(names(logPriors), priornames))

  GeyerThompson <- match.arg(GeyerThompson)
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
    # Circular data in [0, 2pi]
    # NOTE:
    # x <- runif(1000, 0, 2*pi)
    # plot(x, atan2(sin(x), cos(x)))
    # require(circular)
    # plot(x, atan2(sin(circular(x)), cos(circular(x)))) # not enough?
    # z <- atan2(sin(x), cos(x))
    # z[z < 0] <- 2*pi + z[z < 0] # atan2 is right continuous...
    # plot(x, z)
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

  if(phi.p < 1) phi.p <- 1

  ##############################################################
  ######################## Step 2: MCMC ########################
  ##############################################################

  theta.fixo <- log(c(beta.p, phi.p))
  suff.mcmc <- switch(GeyerThompson,
                      large = sampleZ1(beta.p, phi.p, gamma.p, sigma.p, kappa.p,
                                       win, J*100, R_centers, R_clusters,
                                       debug = TRUE),
                      resample = resampleZ1(beta.p, phi.p, final,
                                            win, J*100, R_centers, R_clusters,
                                            debug = TRUE),
                      simple = sampleZ1s(beta.p, phi.p, gamma.p, sigma.p, kappa.p,
                                         win, J, R_centers, R_clusters,
                                         debug = TRUE),
                      importance = sampleZ1(beta.p, phi.p, gamma.p, sigma.p, kappa.p,
                                            win, J, R_centers, R_clusters,
                                            debug = TRUE),
                      MollerAlt = matrix(0, ncol = 2, nrow = J))
  if(verbose) cat(paste0("l = ",1,"\n"))

  # sampled parameters initialization
  beta <- rep(0, L)
  phi <- rep(0, L)
  gamma <- rep(0, L)
  sigma <- rep(0, L)
  kappa <- rep(0, L)

  beta[1] <- beta.p
  phi[1] <- phi.p
  gamma[1] <- gamma.p
  sigma[1] <- sigma.p
  kappa[1] <- kappa.p

  fmat.cond <- sampleThin(fmat, gamma[1], kappa[1], sigma[1], win, R_clusters, up.0)

  # ff <- function(x) posterioriBeta(fmat.cond, beta = x,
  #                                 phi = phi.p, gamma = gamma.p,
  #                                 sigma = sigma.p, kappa = kappa.p,
  #                                 BETA.P = beta.p,
  #                                 suff.beta.phi, theta.fixo, suff.mcmc)
  # # curve(ff, 5, 50)
  # zz <- seq(10, 10000, 1)
  # yy <- numeric(length(zz))
  # for(i in 1:length(zz)) yy[i] <- ff(zz[i])
  # plot(zz,yy)
  suff.all <- sufficientStatAll(fmat.cond, sigma.p, kappa.p,
                                A1 = a1, A2 = a2, B1 = b1, B2 = b2,
                                NN = N, Rc = R_centers, R = R_clusters)
  if(GeyerThompson != "MollerAlt"){
    prob0beta <- posterioriBeta(fmat.cond, beta = beta.p,
                                phi = phi.p, gamma = gamma.p,
                                sigma = sigma.p, kappa = kappa.p,
                                BETA.P = beta.p,
                                suff.beta.phi, theta.fixo, suff.mcmc,
                                SFALL = suff.all,
                                logpriorBeta = logPriors$priorBeta)
    prob0phi <- posterioriPhi(fmat.cond, beta = beta.p,
                              phi = phi.p, gamma = gamma.p,
                              sigma = sigma.p, kapp = kappa.p,
                              PHI.P = phi.p,
                              suff.beta.phi, theta.fixo, suff.mcmc,
                              SFALL = suff.all,
                              logpriorPhi = logPriors$priorPhi)
  } else {
    xAuxiliary0 <- rmat3(tildeTheta[1], tildeTheta[2],
                         gamma.p, sigma.p, kappa.p,
                         win, R_clusters, R_centers)
    suff.beta.phi.x0 <- sufficientStat(xAuxiliary0, Rc = R_centers) # suff.beta.phi
  }
  prob0gamma <- posterioriGamma(fmat.cond, beta = beta.p,
                                phi = phi.p, gamma = gamma.p,
                                sigma = sigma.p, kappa = kappa.p,
                                GAMMA.P = gamma.p,
                                # A1 = a1, A2 = a2, B1 = b1, B2 = b2,
                                # NN = N, Rc = R_centers, R = R_clusters,
                                SFALL = suff.all,
                                logpriorGamma = logPriors$priorGamma)
  prob0sigma <- posterioriSigma(fmat.cond, beta = beta.p,
                                phi = phi.p, gamma = gamma.p,
                                sigma = sigma.p, kappa = kappa.p,
                                SIGMA.P = sigma.p,
                                A1 = a1, A2 = a2, B1 = b1, B2 = b2,
                                NN = N, Rc = R_centers, R = R_clusters,
                                logpriorSigma = logPriors$priorSigma)
  prob0kappa <- posterioriKappa(fmat.cond, beta = beta.p,
                                phi = phi.p, gamma = gamma.p,
                                sigma = sigma.p, kappa = kappa.p,
                                KAPPA.P = kappa.p,
                                A1 = a1, A2 = a2, B1 = b1, B2 = b2,
                                NN = N, Rc = R_centers, R = R_clusters,
                                logpriorKappa = logPriors$priorKappa)

  if(!is.null(resultsName)){
    write.table(t(c(beta[1], phi[1], gamma[1], sigma[1], kappa[1])),
                file = resultsName,
                append = TRUE, row.names = FALSE, col.names = FALSE)
  }

  ini <- 2

  for (l in ini:L) {

    if (l %% 100 == 1) {
      bur <- ifelse(l == 101, 1, 100)
      censorBeta <- min(median(beta[bur:(l-1)]), beta.p)
      censorPhi <- min(median(phi[bur:(l-1)]), 5)
      censorPhi <- max(1/5, censorPhi)
      theta.fixo <- log(c(censorBeta, censorPhi))
      suff.mcmc <- switch(GeyerThompson,
                          large = suff.mcmc, # Mo change
                          resample = suff.mcmc, # Mo change
                          simple = sampleZ1s(censorBeta,
                                             censorPhi,
                                             median(gamma[bur:(l-1)]),
                                             median(sigma[bur:(l-1)]),
                                             median(kappa[bur:(l-1)]),
                                             win, J, R_centers, R_clusters,
                                             debug = TRUE),
                          importance = sampleZ1(censorBeta,
                                                censorPhi,
                                                median(gamma[bur:(l-1)]),
                                                median(sigma[bur:(l-1)]),
                                                median(kappa[bur:(l-1)]),
                                                win, J, R_centers, R_clusters,
                                                debug = TRUE),
                          MollerAlt = suff.mcmc)
    }

    suff.all <- sufficientStatAll(fmat.cond, sigma[l-1], kappa[l-1],
                                  A1 = a1, A2 = a2, B1 = b1, B2 = b2,
                                  NN = N, Rc = R_centers, R = R_clusters)

    if(verbose) cat(paste0("l = ",l,"\n"))

    fmat.cond <- resampleTimes(fmat.cond, R = R_clusters)

    if(GeyerThompson != "MollerAlt"){
      beta.tent <- exp(rnorm(1, log(beta[l-1]), candidateVar[1]))
      prob1 <- posterioriBeta(fmat.cond, beta.tent, phi[l-1], gamma[l-1], sigma[l-1], kappa[l-1],
                              BETA.P = beta.p, suff.beta.phi, theta.fixo, suff.mcmc,
                              SFALL = suff.all,
                              logpriorBeta = logPriors$priorBeta)
      prob0beta <- posterioriBeta(fmat.cond, beta[l-1], phi[l-1], gamma[l-1], sigma[l-1], kappa[l-1],
                                  BETA.P = beta.p, suff.beta.phi, theta.fixo, suff.mcmc,
                                  SFALL = suff.all,
                                  logpriorBeta = logPriors$priorBeta)

      const <- beta.tent/beta[l-1]
      accept1 <- const*exp(prob1-prob0beta)

      if(verbose) cat(paste0("\t ratio beta = ", const, "\t"))
      if(verbose) cat(paste0("\t accept = ", accept1, "\n"))

      if (!is.nan(accept1) && runif(1) < min(1, accept1)) {
        beta[l] <- beta.tent
      } else {
        beta[l] <- beta[l-1]
      }

      # https://stats.stackexchange.com/questions/178226/do-sampling-methods-mcmc-smc-work-for-combination-of-continuous-and-discrete-r
      p.mass <- 0.2 # pnorm(1, log(phi[l-1]), candidateVar[2])
      if(runif(1) < p.mass){
        phi.tent <- 1
      } else {
        phi.tent <- exp(rnorm(1, log(phi[l-1]), candidateVar[2]))
        while(phi.tent < 1) {
          phi.tent <- exp(rnorm(1, log(phi[l-1]), candidateVar[2]))
        }
      }
      prob1 <- posterioriPhi(fmat.cond, beta[l], phi.tent, gamma[l-1], sigma[l-1], kappa[l-1],
                             PHI.P = phi.p,
                             suff.beta.phi, theta.fixo, suff.mcmc,
                             SFALL = suff.all,
                             logpriorPhi = logPriors$priorPhi)
      prob0phi <- posterioriPhi(fmat.cond, beta[l], phi[l-1], gamma[l-1], sigma[l-1], kappa[l-1],
                                PHI.P = phi.p,
                                suff.beta.phi, theta.fixo, suff.mcmc,
                                SFALL = suff.all,
                                logpriorPhi = logPriors$priorPhi)

      const <- phi.tent/phi[l-1]
      accept1 <- const*exp(prob1-prob0phi)

      if(verbose) cat(paste0("\t ratio phi = ", const, "\t"))
      if(verbose) cat(paste0("\t accept = ", accept1, "\n"))

      if (!is.nan(accept1) && runif(1) < min(1, accept1)) {
        phi[l] <- phi.tent
      } else {
        phi[l] <- phi[l-1]
      }
    } else {
      beta.tent <- exp(rnorm(1, log(beta[l-1]), candidateVar[1])) # rlnorm
      p.mass <- 0.2 # pnorm(1, log(phi[l-1]), candidateVar[2])
      if(runif(1) < p.mass){
        phi.tent <- 1
      } else {
        phi.tent <- exp(rnorm(1, log(phi[l-1]), candidateVar[2]))
        while(phi.tent < 1) {
          phi.tent <- exp(rnorm(1, log(phi[l-1]), candidateVar[2]))
        }
      }
      xAuxiliary <- rmat3(beta.tent, phi.tent,
                          gamma[l-1], sigma[l-1], kappa[l-1],
                          win, R_clusters, R_centers)
      suff.beta.phi.x <- sufficientStat(xAuxiliary, Rc = R_centers) # suff.beta.phi

      logH <- (suff.beta.phi.x0[1]-suff.beta.phi.x[1])*log(tildeTheta[1]) + (suff.beta.phi.x0[2]-suff.beta.phi.x[2])*log(tildeTheta[2])
      logH <- logH + (suff.beta.phi.x0[1]-suff.beta.phi[1])*log(beta[l-1]) + (suff.beta.phi.x0[2]-suff.beta.phi[2])*log(phi[l-1])
      logH <- logH - (suff.beta.phi.x[1]-suff.beta.phi[1])*log(beta.tent) - (suff.beta.phi.x[2]-suff.beta.phi[2])*log(phi.tent)
      logH <- logH + logPriors$priorBeta(beta.tent, beta.p) - logPriors$priorBeta(beta[l-1], beta.p)
      logH <- logH + logPriors$priorPhi(phi.tent, phi.p) - logPriors$priorPhi(phi[l-1], phi.p)
      logH <- logH + dlnorm(beta[l-1], log(beta.tent), candidateVar[1]) - dlnorm(beta.tent, log(beta[l-1]), candidateVar[1])
      logH <- logH + dlnorm(phi[l-1], log(phi.tent), candidateVar[2]) - dlnorm(phi.tent, log(phi[l-1]), candidateVar[2])

      if (!is.nan(logH) && runif(1) < min(1, exp(logH))) {
        beta[l] <- beta.tent
        phi[l] <- phi.tent
        xAuxiliary0 <- xAuxiliary
      } else {
        beta[l] <- beta[l-1]
        phi[l] <- phi[l-1]
      }
    }

    gamma.tent <- exp(rnorm(1, log(gamma[l-1]), candidateVar[3]))
    prob1 <- posterioriGamma(fmat.cond, beta[l], phi[l], gamma.tent, sigma[l-1], kappa[l-1],
                             GAMMA.P = gamma.p,
                             SFALL = suff.all,
                             logpriorGamma = logPriors$priorGamma)
    prob0gamma <- posterioriGamma(fmat.cond, beta[l], phi[l], gamma[l-1], sigma[l-1], kappa[l-1],
                                  GAMMA.P = gamma.p,
                                  SFALL = suff.all,
                                  logpriorGamma = logPriors$priorGamma)

    const <- gamma.tent/gamma[l-1]
    accept1 <- const*exp(prob1-prob0gamma)

    if(verbose) cat(paste0("\t ratio gamma = ", const, "\t"))
    if(verbose) cat(paste0("\t accept = ", accept1, "\n"))

    if (!is.nan(accept1) && runif(1) < min(1, accept1)) {
      gamma[l] <- gamma.tent
    } else {
      gamma[l] <- gamma[l-1]
    }

    sigma.tent <- exp(rnorm(1, log(sigma[l-1]), candidateVar[4]))
    prob1 <- posterioriSigma(fmat.cond, beta[l], phi[l], gamma[l], sigma.tent, kappa[l-1],
                             SIGMA.P = sigma.p,
                             A1 = a1, A2 = a2, B1 = b1, B2 = b2,
                             NN = N, Rc = R_centers, R = R_clusters,
                             logpriorSigma = logPriors$priorSigma)
    prob0sigma <- posterioriSigma(fmat.cond, beta[l], phi[l], gamma[l], sigma[l-1], kappa[l-1],
                                  SIGMA.P = sigma.p,
                                  A1 = a1, A2 = a2, B1 = b1, B2 = b2,
                                  NN = N, Rc = R_centers, R = R_clusters,
                                  logpriorSigma = logPriors$priorSigma)

    const <- sigma.tent/sigma[l-1]
    accept1 <- const*exp(prob1-prob0sigma)

    if(verbose) cat(paste0("\t ratio sigma = ", const, "\t"))
    if(verbose) cat(paste0("\t accept = ", accept1, "\n"))

    if (!is.nan(accept1) && runif(1) < min(1, accept1)) {
      sigma[l] <- sigma.tent
    } else {
      sigma[l] <- sigma[l-1]
    }

    kappa.tent <- exp(rnorm(1, log(kappa[l-1]), candidateVar[5]))
    prob1 <- posterioriKappa(fmat.cond, beta[l], phi[l], gamma[l], sigma[l], kappa.tent,
                             KAPPA.P = kappa.p,
                             A1 = a1, A2 = a2, B1 = b1, B2 = b2,
                             NN = N, Rc = R_centers, R = R_clusters,
                             logpriorKappa = logPriors$priorKappa)
    prob0kappa <- posterioriKappa(fmat.cond, beta[l], phi[l], gamma[l], sigma[l], kappa[l-1],
                                  KAPPA.P = kappa.p,
                                  A1 = a1, A2 = a2, B1 = b1, B2 = b2,
                                  NN = N, Rc = R_centers, R = R_clusters,
                                  logpriorKappa = logPriors$priorKappa)

    const <- kappa.tent/kappa[l-1]
    accept1 <- const*exp(prob1-prob0kappa)

    if(verbose) cat(paste0("\t ratio kappa = ", const, "\t"))
    if(verbose) cat(paste0("\t accept = ", accept1, "\n"))

    if (!is.nan(accept1) && runif(1) < min(1, accept1)) {
      kappa[l] <- kappa.tent
    } else {
      kappa[l] <- kappa[l-1]
    }

    if(!is.null(resultsName)){
      write.table(t(c(beta[l], phi[l], gamma[l], sigma[l], kappa[l])),
                  file = resultsName,
                  append = TRUE, row.names = FALSE, col.names = FALSE)
    }

    fmat.cond <- sampleThin(fmat.cond, gamma[l], kappa[l], sigma[l],
                            win, R_clusters, up.0)

  }

  parameters <- cbind(beta, phi, gamma, sigma, kappa)
  colnames(parameters) <- c("beta","phi","gamma","sigma","kappa")

  results <- list(parameters = parameters,
                  augmentedData = fmat.cond,
                  window = win,
                  R_centers = R_centers,
                  R_clusters = R_clusters,
                  suffMCMC = suff.mcmc,
                  logPriors = logPriors)
  class(results) <- c("fitmat3", "list")

  return(results)

}
