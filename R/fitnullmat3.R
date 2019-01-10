#' Fit a Poisson marked point process model via MCMC sampling
#'
#' Markov chain Monte Carlo sampler from Poisson marked spatial point process
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
#' @param L Number of samples for the MCMC sampler. Defaults to 10000.
#' @param NN Number of samples to Monte-Carlo estimation of constant. Defaults to 2000.
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
#' @param candidateVar Variability of MCMC candidate selection step, if tuning is necessary.
#' @param verbose Prints acceptance rates for candidates. Defaults to TRUE.
#' @param ... further arguments passed to \code{\link{read.table}}.
#'
#' @export
#' @useDynLib mat3c
#' @import circular spatstat
#' @importFrom Rcpp evalCpp
#' @return \code{fitnullmat3} returns a list containing at least the following components
#'   \item{parameters}{An L by 5 matrix with the samples from beta, phi, gamma, sigma, kappa.}
#'
#' @examples
#' set.seed(1234)
#' x <- rmat3(70, 2, 5, 0.05, 3)
#' plot(x)
#' # Changing default sampling sizes to make it run fast
#' model <- fitnullmat3(x, L = 100, seed = 1234)
#'
#' @author Guilherme Ludwig and Nancy Garcia
#'
#' @references
#'
#'   Garcia, N., Guttorp, P. and Ludwig, G. (2018) TBD
#'
#'   Geyer, C. J., and Thompson, E. A. (1992) "Constrained Monte Carlo maximum likelihood
#'   for dependent data." Journal of the Royal Statistical Society. Series B (Methodological),
#'   657-699.
#'
#' @seealso \code{\link{rmat3}}
#' @keywords Spatial Statistics
#' @keywords Functional Data Analysis
fitnullmat3 <- function(mat3, fname = c("centers.txt", "fingers.txt"),
                    win = owin(c(0,1), c(0,1)),
                    L = 10000, NN = 2000, seed = NULL,
                    resultsName = NULL,
                    logPriors = preparePriors(),
                    initialValues = pickInitialValues(),
                    candidateVar = c(0.1, 0.2, 0.1, 0.1, 0.1),
                    verbose = TRUE,
                    ...){

  set.seed(seed)

  priornames <- c("priorBeta", "priorPhi", "priorGamma", "priorSigma", "priorKappa")
  stopifnot(all.equal(names(logPriors), priornames))

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

  suff.beta.phi <- c(length(final), 0) # sufficientStat(final, Rc = R_centers)

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
  phi.p   <- 1 # Does not exist
  kappa.p <- initialValues$iniKappa(angles)
  sigma.p <- initialValues$iniSigma(fmat)

  if(phi.p < 1) phi.p <- 1

  ##############################################################
  ######################## Step 2: MCMC ########################
  ##############################################################

  theta.fixo <- log(c(beta.p, phi.p))

  if(verbose) cat(paste0("l = ",1,"\n"))

  # sampled parameters initialization
  beta <- rep(0, L)
  phi <- rep(0, L)
  gamma <- rep(0, L)
  sigma <- rep(0, L)
  kappa <- rep(0, L)

  beta[1] <- beta.p
  phi[1] <- 1
  gamma[1] <- gamma.p
  sigma[1] <- sigma.p
  kappa[1] <- kappa.p

  # fmat.cond <- sampleThin(fmat, gamma[1], kappa[1], sigma[1], win, R_clusters, up.0)

  z3 <- sum(sapply(final, function(x) nrow(x$fingers)))
  z5 <- 0
  for(i in 1:length(final)){
    z1 <- final[[i]]$fingers[,1] - final[[i]]$centers[1]
    z2 <- final[[i]]$fingers[,2] - final[[i]]$centers[2]
    radius <- z1^2 + z2^2
    angle <- atan2(z2, z1)
    angle[angle < 0] <- 2*pi + angle[angle < 0] # atan2 -> [-pi,pi] but we sample u on [0, 2pi]
    direction <- mean(circular(angle, modulo = "2pi"))
    X.st <- rnorm(NN, 0, 1)
    temp1 <- circular:::RvonmisesRad(NN, direction, kappa[1])
    X <- cbind(final[[i]]$centers[1] + sigma[1]*abs(X.st)*cos(temp1),
               final[[i]]$centers[2] + sigma[1]*abs(X.st)*sin(temp1))
    d.s.1 <- !(X[, 1] > a1 & X[, 1] < a2 & X[, 2] > b1 & X[, 2] < b2)
    n.s.1 <- sum(d.s.1)
    prob.tot <- 1 - n.s.1/NN
    z5 <- z5 + sum(circular:::DvonmisesRad(angle, direction, kappa[1], log = TRUE)
                   + 0.5*(log(2/pi) + log(radius)) - log(sigma[1]) - 0.5*radius/(sigma[1]^2) -
                     log(prob.tot))
    # z5 <- sum(outer(circular:::DvonmisesRad(seq(0, 2*pi - pi/20, pi/20), runif(1, 0, 2*pi), kappa[1]),
    #                 dnorm(seq(-3,3,.1)*sigma[1], 0, sigma[1])))*(pi/20)*.1*sigma[1]
  }

  suff.all <- c(length(final), 0, z3, 0, z5)
  # suff.all <- sufficientStatAll(fmat.cond, sigma.p, kappa.p,
  #                               A1 = a1, A2 = a2, B1 = b1, B2 = b2,
  #                               NN = N, Rc = R_centers, R = R_clusters)

  # posterioriBeta
  theta.prime <- log(c(beta[1], phi[1])) # Candidates
  prob0beta <- sum(suff.beta.phi*theta.prime) - gamma[1]*suff.all[1] + log(gamma[1])*suff.all[3] -
    suff.all[1]*log(1 - exp(-gamma[1])) + logPriors$priorBeta(beta[1], beta.p)

  # posterioriPhi
  prob0phi <- 1

  # posterioriGamma
  prob0gamma <- -gamma[1]*suff.all[1] + log(gamma[1])*suff.all[3] + suff.all[5] - suff.all*log(1 - exp(-gamma[1]))
  prob0gamma <- prob0gamma + logPriors$priorGamma(gamma[1], gamma.p)

  # posterioriSigma
  z5 <- 0
  for(i in 1:length(final)){
    z1 <- final[[i]]$fingers[,1] - final[[i]]$centers[1]
    z2 <- final[[i]]$fingers[,2] - final[[i]]$centers[2]
    radius <- z1^2 + z2^2
    angle <- atan2(z2, z1)
    angle[angle < 0] <- 2*pi + angle[angle < 0]
    direction <- mean(circular(angle, modulo = "2pi"))
    X.st <- rnorm(NN, 0, 1)
    temp1 <- circular:::RvonmisesRad(NN, direction, kappa[1])
    X <- cbind(final[[i]]$centers[1] + sigma[1]*abs(X.st)*cos(temp1),
               final[[i]]$centers[2] + sigma[1]*abs(X.st)*sin(temp1))
    d.s.1 <- !(X[, 1] > a1 & X[, 1] < a2 & X[, 2] > b1 & X[, 2] < b2)
    n.s.1 <- sum(d.s.1)
    prob.tot <- 1 - n.s.1/NN
    z5 <- z5 + sum(circular:::DvonmisesRad(angle, direction, kappa[1], log = TRUE)
                   + 0.5*(log(2/pi) + log(radius)) - log(sigma[1]) - 0.5*radius/(sigma[1]^2) -
                     log(prob.tot))
  }
  suff.all <- c(length(final), 0, z3, 0, z5)
  prob0sigma <- -gamma[1]*suff.all[1] + log(gamma[1])*suff.all[3] + suff.all[5] - suff.all[1]*log(1 - exp(-gamma[1]))
  prob0sigma <- prob0sigma + logPriors$priorSigma(sigma[1], sigma.p)

  # posterioriKappa
  z5 <- 0
  for(i in 1:length(final)){
    z1 <- final[[i]]$fingers[,1] - final[[i]]$centers[1]
    z2 <- final[[i]]$fingers[,2] - final[[i]]$centers[2]
    radius <- z1^2 + z2^2
    angle <- atan2(z2, z1)
    angle[angle < 0] <- 2*pi + angle[angle < 0]
    direction <- mean(circular(angle, modulo = "2pi"))
    X.st <- rnorm(NN, 0, 1)
    temp1 <- circular:::RvonmisesRad(NN, direction, kappa[1])
    X <- cbind(final[[i]]$centers[1] + sigma[1]*abs(X.st)*cos(temp1),
               final[[i]]$centers[2] + sigma[1]*abs(X.st)*sin(temp1))
    d.s.1 <- !(X[, 1] > a1 & X[, 1] < a2 & X[, 2] > b1 & X[, 2] < b2)
    n.s.1 <- sum(d.s.1)
    prob.tot <- 1 - n.s.1/NN
    z5 <- z5 + sum(circular:::DvonmisesRad(angle, direction, kappa[1], log = TRUE)
                   + 0.5*(log(2/pi) + log(radius)) - log(sigma[1]) - 0.5*radius/(sigma[1]^2) -
                     log(prob.tot))
  }
  suff.all <- c(length(final), 0, z3, 0, z5)
  prob0kappa <- -gamma[1]*suff.all[1] + log(gamma[1])*suff.all[3] + suff.all[5] - suff.all[1]*log(1 - exp(-gamma[1]))
  prob0kappa <- prob0kappa + logPriors$priorKappa(kappa[1], kappa.p)

  if(!is.null(resultsName)){
    write.table(t(c(beta[1], phi[1], gamma[1], sigma[1], kappa[1])),
                file = resultsName,
                append = TRUE, row.names = FALSE, col.names = FALSE)
  }

  ini <- 2

  for (l in ini:L) {

    z3 <- sum(sapply(final, function(x) nrow(x$fingers)))
    z5 <- 0
    for(i in 1:length(final)){
      z1 <- final[[i]]$fingers[,1] - final[[i]]$centers[1]
      z2 <- final[[i]]$fingers[,2] - final[[i]]$centers[2]
      radius <- z1^2 + z2^2
      angle <- atan2(z2, z1)
      angle[angle < 0] <- 2*pi + angle[angle < 0]
      direction <- mean(circular(angle, modulo = "2pi"))
      X.st <- rnorm(NN, 0, 1)
      temp1 <- circular:::RvonmisesRad(NN, direction, kappa[l-1])
      X <- cbind(final[[i]]$centers[1] + sigma[l-1]*abs(X.st)*cos(temp1),
                 final[[i]]$centers[2] + sigma[l-1]*abs(X.st)*sin(temp1))
      d.s.1 <- !(X[, 1] > a1 & X[, 1] < a2 & X[, 2] > b1 & X[, 2] < b2)
      n.s.1 <- sum(d.s.1)
      prob.tot <- 1 - n.s.1/NN
      z5 <- z5 + sum(circular:::DvonmisesRad(angle, direction, kappa[l-1], log = TRUE)
                     + 0.5*(log(2/pi) + log(radius)) - log(sigma[l-1]) - 0.5*radius/(sigma[l-1]^2) -
                       log(prob.tot))
    }
    suff.all <- c(length(final), 0, z3, 0, z5)

    if(verbose) cat(paste0("l = ",l,"\n"))

    # fmat.cond <- resampleTimes(fmat.cond, R = R_clusters)

    beta.tent <- exp(rnorm(1, log(beta[l-1]), candidateVar[1]))
    # prob1
    theta.prime <- log(c(beta.tent, phi[l-1])) # Candidates
    prob1 <- sum(suff.beta.phi*theta.prime) - gamma[l-1]*suff.all[1] + log(gamma[l-1])*suff.all[3] -
      suff.all[1]*log(1 - exp(-gamma[l-1])) + logPriors$priorBeta(beta.tent, beta.p)
    # prob0beta
    theta.prime <- log(c(beta[l-1], phi[l-1])) # Candidates
    prob0beta <- sum(theta.prime*suff.beta.phi) - gamma[l-1]*suff.all[1] + log(gamma[l-1])*suff.all[3] -
      suff.all[1]*log(1 - exp(-gamma[l-1])) + logPriors$priorBeta(beta[l-1], beta.p)

    const <- beta.tent/beta[l-1]
    accept1 <- const*exp(prob1-prob0beta)

    if(verbose) cat(paste0("\t ratio beta = ", const, "\t"))
    if(verbose) cat(paste0("\t accept = ", accept1, "\n"))

    if (!is.nan(accept1) && runif(1) < min(1, accept1)) {
      beta[l] <- beta.tent
    } else {
      beta[l] <- beta[l-1]
    }

    # phi
    phi[l] <- 1

    # gamma
    gamma.tent <- exp(rnorm(1, log(gamma[l-1]), candidateVar[3]))

    prob1 <- -gamma.tent*suff.all[1] + log(gamma.tent)*suff.all[3] + suff.all[5] - suff.all*log(1 - exp(-gamma.tent))
    prob1 <- prob1 + logPriors$priorGamma(gamma.tent, gamma.p)
    prob0gamma <- -gamma[l-1]*suff.all[1] + log(gamma[l-1])*suff.all[3] + suff.all[5] - suff.all*log(1 - exp(-gamma[l-1]))
    prob0gamma <- prob0gamma + logPriors$priorGamma(gamma[l-1], gamma.p)

    const <- gamma.tent/gamma[l-1]
    accept1 <- const*exp(prob1-prob0gamma)

    if(verbose) cat(paste0("\t ratio gamma = ", const, "\t"))
    if(verbose) cat(paste0("\t accept = ", accept1, "\n"))

    if (!is.nan(accept1) && runif(1) < min(1, accept1)) {
      gamma[l] <- gamma.tent
    } else {
      gamma[l] <- gamma[l-1]
    }

    # sigma
    sigma.tent <- exp(rnorm(1, log(sigma[l-1]), candidateVar[4]))
    z5 <- 0
    for(i in 1:length(final)){
      z1 <- final[[i]]$fingers[,1] - final[[i]]$centers[1]
      z2 <- final[[i]]$fingers[,2] - final[[i]]$centers[2]
      radius <- z1^2 + z2^2
      angle <- atan2(z2, z1)
      angle[angle < 0] <- 2*pi + angle[angle < 0]
      direction <- mean(circular(angle, modulo = "2pi"))
      X.st <- rnorm(NN, 0, 1)
      temp1 <- circular:::RvonmisesRad(NN, direction, kappa[l-1])
      X <- cbind(final[[i]]$centers[1] + sigma.tent*abs(X.st)*cos(temp1),
                 final[[i]]$centers[2] + sigma.tent*abs(X.st)*sin(temp1))
      d.s.1 <- !(X[, 1] > a1 & X[, 1] < a2 & X[, 2] > b1 & X[, 2] < b2)
      n.s.1 <- sum(d.s.1)
      prob.tot <- 1 - n.s.1/NN
      z5 <- z5 + sum(circular:::DvonmisesRad(angle, direction, kappa[l-1], log = TRUE)
                     + 0.5*(log(2/pi) + log(radius)) - log(sigma.tent) - 0.5*radius/(sigma.tent^2) -
                       log(prob.tot))
    }
    suff.all <- c(length(final), 0, z3, 0, z5)
    prob1 <- -gamma[l]*suff.all[1] + log(gamma[l])*suff.all[3] + suff.all[5] - suff.all[1]*log(1 - exp(-gamma[l]))
    prob1 <- prob1 + logPriors$priorSigma(sigma.tent, sigma.p)
    z5 <- 0
    for(i in 1:length(final)){
      z1 <- final[[i]]$fingers[,1] - final[[i]]$centers[1]
      z2 <- final[[i]]$fingers[,2] - final[[i]]$centers[2]
      radius <- z1^2 + z2^2
      angle <- atan2(z2, z1)
      angle[angle < 0] <- 2*pi + angle[angle < 0]
      direction <- mean(circular(angle, modulo = "2pi"))
      X.st <- rnorm(NN, 0, 1)
      temp1 <- circular:::RvonmisesRad(NN, direction, kappa[l-1])
      X <- cbind(final[[i]]$centers[1] + sigma[l-1]*abs(X.st)*cos(temp1),
                 final[[i]]$centers[2] + sigma[l-1]*abs(X.st)*sin(temp1))
      d.s.1 <- !(X[, 1] > a1 & X[, 1] < a2 & X[, 2] > b1 & X[, 2] < b2)
      n.s.1 <- sum(d.s.1)
      prob.tot <- 1 - n.s.1/NN
      z5 <- z5 + sum(circular:::DvonmisesRad(angle, direction, kappa[l-1], log = TRUE)
                     + 0.5*(log(2/pi) + log(radius)) - log(sigma[l-1]) - 0.5*radius/(sigma[l-1]^2) -
                       log(prob.tot))
    }
    suff.all <- c(length(final), 0, z3, 0, z5)
    prob0sigma <- -gamma[l]*suff.all[1] + log(gamma[l])*suff.all[3] + suff.all[5] - suff.all[1]*log(1 - exp(-gamma[l]))
    prob0sigma <- prob0sigma + logPriors$priorSigma(sigma[l-1], sigma.p)

    const <- sigma.tent/sigma[l-1]
    accept1 <- const*exp(prob1-prob0sigma)

    if(verbose) cat(paste0("\t ratio sigma = ", const, "\t"))
    if(verbose) cat(paste0("\t accept = ", accept1, "\n"))

    if (!is.nan(accept1) && runif(1) < min(1, accept1)) {
      sigma[l] <- sigma.tent
    } else {
      sigma[l] <- sigma[l-1]
    }

    # posterioriKappa
    kappa.tent <- exp(rnorm(1, log(kappa[l-1]), candidateVar[5]))
    z5 <- 0
    for(i in 1:length(final)){
      z1 <- final[[i]]$fingers[,1] - final[[i]]$centers[1]
      z2 <- final[[i]]$fingers[,2] - final[[i]]$centers[2]
      radius <- z1^2 + z2^2
      angle <- atan2(z2, z1)
      angle[angle < 0] <- 2*pi + angle[angle < 0]
      direction <- mean(circular(angle, modulo = "2pi"))
      X.st <- rnorm(NN, 0, 1)
      temp1 <- circular:::RvonmisesRad(NN, direction, kappa.tent)
      X <- cbind(final[[i]]$centers[1] + sigma[l]*abs(X.st)*cos(temp1),
                 final[[i]]$centers[2] + sigma[l]*abs(X.st)*sin(temp1))
      d.s.1 <- !(X[, 1] > a1 & X[, 1] < a2 & X[, 2] > b1 & X[, 2] < b2)
      n.s.1 <- sum(d.s.1)
      prob.tot <- 1 - n.s.1/NN
      z5 <- z5 + sum(circular:::DvonmisesRad(angle, direction, kappa.tent, log = TRUE)
                     + 0.5*(log(2/pi) + log(radius)) - log(sigma[l]) - 0.5*radius/(sigma[l]^2) -
                       log(prob.tot))
    }
    suff.all <- c(length(final), 0, z3, 0, z5)
    prob0kappa <- -gamma[l]*suff.all[1] + log(gamma[l])*suff.all[3] + suff.all[5] - suff.all[1]*log(1 - exp(-gamma[l]))
    prob0kappa <- prob0kappa + logPriors$priorKappa(kappa[l-1], kappa.p)

    z5 <- 0
    for(i in 1:length(final)){
      z1 <- final[[i]]$fingers[,1] - final[[i]]$centers[1]
      z2 <- final[[i]]$fingers[,2] - final[[i]]$centers[2]
      radius <- z1^2 + z2^2
      angle <- atan2(z2, z1)
      angle[angle < 0] <- 2*pi + angle[angle < 0]
      direction <- mean(circular(angle, modulo = "2pi"))
      X.st <- rnorm(NN, 0, 1)
      temp1 <- circular:::RvonmisesRad(NN, direction, kappa[l-1])
      X <- cbind(final[[i]]$centers[1] + sigma[l]*abs(X.st)*cos(temp1),
                 final[[i]]$centers[2] + sigma[l]*abs(X.st)*sin(temp1))
      d.s.1 <- !(X[, 1] > a1 & X[, 1] < a2 & X[, 2] > b1 & X[, 2] < b2)
      n.s.1 <- sum(d.s.1)
      prob.tot <- 1 - n.s.1/NN
      z5 <- z5 + sum(circular:::DvonmisesRad(angle, direction, kappa[l-1], log = TRUE)
                     + 0.5*(log(2/pi) + log(radius)) - log(sigma[l]) - 0.5*radius/(sigma[l]^2) -
                       log(prob.tot))
    }
    suff.all <- c(length(final), 0, z3, 0, z5)
    prob1 <- -gamma[l]*suff.all[1] + log(gamma[l])*suff.all[3] + suff.all[5] - suff.all[1]*log(1 - exp(-gamma[l]))
    prob1 <- prob1 + logPriors$priorKappa(kappa.tent, kappa.p)

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

    # fmat.cond <- sampleThin(fmat.cond, gamma[l], kappa[l], sigma[l],
    #                         win, R_clusters, up.0)

  }

  parameters <- cbind(beta, phi, gamma, sigma, kappa)
  colnames(parameters) <- c("beta","phi","gamma","sigma","kappa")

  results <- list(parameters = parameters,
                  # augmentedData = fmat.cond,
                  window = win,
                  # R_centers = R_centers,
                  # R_clusters = R_clusters,
                  # suffMCMC = suff.mcmc,
                  logPriors = logPriors)
  class(results) <- c("fitmat3", "list")

  return(results)

}
