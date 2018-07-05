plotPriors <- function(object){
  win <- object$window
  a1 <- win$x[1]
  a2 <- win$x[2]
  b1 <- win$y[1]
  b2 <- win$y[2]
  limx <- c(a1, a2)
  limy <- c(b1, b2)

  stop() # TODO

  theta.fixo <- log(c(object$parameters[1,1:2]))
  suff.mcmc <- object$suffMCMC
  beta.p <- object$parameters[1,1]
  phi.p <- object$parameters[1,2]
  gamma.p <- object$parameters[1,3]
  sigma.p <- object$parameters[1,4]
  kappa.p <- object$parameters[1,5]

  R_centers <- object$R_centers
  R_clusters <- object$R_clusters

  fmat.cond <- object$augmentedData
  suff.beta.phi <- sufficientStat(final, R = R_centers)

  f_beta <- function(x) mat3c:::posterioriBeta(fmat.cond, beta = x,
                                               phi = phi.p, gamma = gamma.p,
                                               sigma = sigma.p, kappa = kappa.p,
                                               BETA.P = beta.p,
                                               suff.beta.phi, theta.fixo, suff.mcmc)
  f_phi <- function(x) mat3c:::posterioriPhi(fmat.cond, beta = beta.p,
                                             phi = x, gamma = gamma.p,
                                             sigma = sigma.p, kapp = kappa.p,
                                             PHI.P = phi.p,
                                             suff.beta.phi, theta.fixo, suff.mcmc)

  fmat.cond <- fmat
  N <- 10000
  f_gamma <- function(x) mat3c:::posterioriGamma(fmat.cond, beta = beta.p,
                                                 phi = phi.p, gamma = x,
                                                 sigma = sigma.p, kappa = kappa.p,
                                                 GAMMA.P = gamma.p,
                                                 A1 = a1, A2 = a2, B1 = b1, B2 = b2,
                                                 NN = N, Rc = R_centers, R = R_clusters)
  f_sigma <- function(x) mat3c:::posterioriSigma(fmat.cond, beta = beta.p,
                                                 phi = phi.p, gamma = gamma.p,
                                                 sigma = x, kappa = kappa.p,
                                                 SIGMA.P = sigma.p,
                                                 A1 = a1, A2 = a2, B1 = b1, B2 = b2,
                                                 NN = N, Rc = R_centers, R = R_clusters)
  f_kappa <- function(x) mat3c:::posterioriKappa(fmat.cond, beta = beta.p,
                                                 phi = phi.p, gamma = gamma.p,
                                                 sigma = sigma.p, kappa = x,
                                                 KAPPA.P = kappa.p,
                                                 A1 = a1, A2 = a2, B1 = b1, B2 = b2,
                                                 NN = N, Rc = R_centers, R = R_clusters)

  I <- 200
  beta0 <- seq(5, 250, length.out = I) # 125
  phi0 <- seq(0.5, 5, length.out = I) # 2
  gamma0 <- seq(10, 100,length.out = I) # 50
  sigma0 <- seq(0.01, 0.05, length.out = I) # 0.025
  kappa0 <- seq(0.5, 6, length.out = I) # 3

  betaPos <- betaPrior <- betaLik <- numeric(I)
  phiPos <- phiPrior <- phiLik <- numeric(I)
  gammaPos <- gammaPrior <- gammaLik <- numeric(I)
  sigmaPos <- sigmaPrior <- sigmaLik <- numeric(I)
  kappaPos <- kappaPrior <- kappaLik <- numeric(I)

  for(i in seq_len(I)){
    print(i)
    # BETA
    betaPos[i] <- f_beta(beta0[i])
    betaPrior[i] <- dlnorm(beta0[i], log(beta.p), 1, log = TRUE)
    betaLik[i] <- betaPos[i] - betaPrior[i]
    # PHI
    phiPos[i] <- f_phi(phi0[i])
    phiPrior[i] <- dgamma(phi0[i], shape = 8/3, scale = 3/4, log = TRUE)
    phiLik[i] <- phiPos[i] - phiPrior[i]
    # GAMMA
    gammaPos[i] <- f_gamma(gamma0[i])
    gammaPrior[i] <- dgamma(gamma0[i], shape = gamma.p, rate = 1, log = TRUE)
    gammaLik[i] <- gammaPos[i] - gammaPrior[i]
    # SIGMA
    sigmaPos[i] <- f_sigma(sigma0[i])
    alpha <- (sigma.p^2)/10
    sigmaPrior[i] <- alpha * log(sigma.p) - lgamma(alpha) - (alpha + 1) * log(sigma0[i]) - (sigma.p/sigma0[i])
    sigmaLik[i] <- sigmaPos[i] - sigmaPrior[i]
    # KAPPA
    kappaPos[i] <- f_kappa(kappa0[i])
    kappaPrior[i] <- dgamma(kappa0[i], shape = (kappa.p^2)/10, scale = 10/kappa.p, log = TRUE)
    kappaLik[i] <- kappaPos[i] - kappaPrior[i]
  }

  layout(matrix(c(1:6), ncol=2))
  # BETA
  plot(beta0, betaPos, type = "l", ylim = range(c(betaPos,betaPrior,betaLik)))
  lines(beta0, betaLik, col="Blue")
  abline(v = 125, col="Blue", lty = 2)
  lines(beta0, betaPrior, col="Red")
  abline(v = beta.p, col="Red", lty = 2)
  # PHI
  plot(phi0, phiPos, type = "l", ylim = range(c(phiPos,phiPrior,phiLik)))
  lines(phi0, phiLik, col="Blue")
  abline(v = 2, col="Blue", lty = 2)
  lines(phi0, phiPrior, col="Red")
  abline(v = phi.p, col="Red", lty = 2)
  # GAMMA
  plot(gamma0, gammaPos, type = "l", ylim = range(c(gammaPos,gammaPrior,gammaLik)))
  lines(gamma0, gammaLik, col="Blue")
  abline(v = 50, col="Blue", lty = 2)
  lines(gamma0, gammaPrior, col="Red")
  abline(v = gamma.p, col="Red", lty = 2)
  # SIGMA
  plot(sigma0, sigmaPos, type = "l", ylim = range(c(sigmaPos,sigmaPrior,sigmaLik)))
  lines(sigma0, sigmaLik, col="Blue")
  abline(v = 0.025, col="Blue", lty = 2)
  lines(sigma0, sigmaPrior, col="Red")
  abline(v = sigma.p, col="Red", lty = 2)
  # KAPPA
  plot(kappa0, kappaPos, type = "l", ylim = range(c(kappaPos,kappaPrior,kappaLik)))
  lines(kappa0, kappaLik, col="Blue")
  abline(v = 3, col="Blue", lty = 2)
  lines(kappa0, kappaPrior, col="Red")
  abline(v = kappa.p, col="Red", lty = 2)
  # LEGEND
  plot(0, type="n", axes=FALSE, xlab="", ylab="")
  legend("center", lty=c(1,1,1,2,2), col=c("Black", "Blue", "Red","Blue", "Red"),
         legend = c("Posterior", "Likelihood", "Prior", "Truth", "Starting value"))
  ############
  layout(matrix(c(1:6), ncol=2))
  # BETA
  plot(beta0, exp(1)^betaPos, type = "l", ylim = range(exp(1)^c(betaPos,betaPrior,betaLik)))
  lines(beta0, exp(1)^betaLik, col="Blue")
  abline(v = 125, col="Blue", lty = 2)
  lines(beta0, exp(1)^betaPrior, col="Red")
  abline(v = beta.p, col="Red", lty = 2)
  # PHI
  plot(phi0, exp(1)^phiPos, type = "l", ylim = range(exp(1)^c(phiPos,phiPrior,phiLik)))
  lines(phi0, exp(1)^phiLik, col="Blue")
  abline(v = 2, col="Blue", lty = 2)
  lines(phi0, exp(1)^phiPrior, col="Red")
  abline(v = phi.p, col="Red", lty = 2)
  # GAMMA
  plot(gamma0, exp(1)^gammaPos, type = "l", ylim = range(exp(1)^c(gammaPos,gammaPrior,gammaLik)))
  lines(gamma0, exp(1)^gammaLik, col="Blue")
  abline(v = 50, col="Blue", lty = 2)
  lines(gamma0, exp(1)^gammaPrior, col="Red")
  abline(v = gamma.p, col="Red", lty = 2)
  # SIGMA
  plot(sigma0, exp(1)^sigmaPos, type = "l", ylim = range(exp(1)^c(sigmaPos,sigmaPrior,sigmaLik)))
  lines(sigma0, exp(1)^sigmaLik, col="Blue")
  abline(v = 0.025, col="Blue", lty = 2)
  lines(sigma0, exp(1)^sigmaPrior, col="Red")
  abline(v = sigma.p, col="Red", lty = 2)
  # KAPPA
  plot(kappa0, exp(1)^kappaPos, type = "l", ylim = range(exp(1)^c(kappaPos,kappaPrior,kappaLik)))
  lines(kappa0, exp(1)^kappaLik, col="Blue")
  abline(v = 3, col="Blue", lty = 2)
  lines(kappa0, exp(1)^kappaPrior, col="Red")
  abline(v = kappa.p, col="Red", lty = 2)
  # LEGEND
  plot(0, type="n", axes=FALSE, xlab="", ylab="")
  legend("center", lty=c(1,1,1,2,2), col=c("Black", "Blue", "Red","Blue", "Red"),
         legend = c("Posterior", "Likelihood", "Prior", "Truth", "Starting value"))
}
