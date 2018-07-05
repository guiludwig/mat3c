posterioriBeta <- function(proc.fmat, beta, phi, gamma, sigma, kappa,
                           BETA.P = beta.p,
                           SUFF.BETA.PHI = suff.beta.phi,
                           THETA.FIXO = theta.fixo, # holds parameters from which suff.mcmc is sampled
                           SUFF.MCMC = suff.mcmc,
                           SFALL = suff.all,
                           logpriorBeta){
  theta.prime <- log(c(beta, phi))
  part1 <- sum(theta.prime*SUFF.BETA.PHI)
  JJ <- nrow(SUFF.MCMC)
  temp <- numeric(JJ)
  for (j in 1:JJ){
    temp[j] <- sum((theta.prime-THETA.FIXO)*SUFF.MCMC[j, ]) # < theta - psi, T*(N) >
  }
  part2 <- mean(exp(temp))
  # like.theta <- (part1 - log(part2))
  like.theta <- (part1 - log(part2)) - gamma*SFALL[1] + log(gamma)*SFALL[3] + gamma*SFALL[4] -
    SFALL[1]*log(1 - exp(-gamma))
  return(logpriorBeta(beta, BETA.P) + like.theta)
}

posterioriPhi <- function(proc.fmat, beta, phi, gamma, sigma, kappa,
                          PHI.P = phi.p,
                          SUFF.BETA.PHI = suff.beta.phi,
                          THETA.FIXO = theta.fixo, # holds beta.p, phi.p
                          SUFF.MCMC = suff.mcmc,
                          SFALL = suff.all,
                          logpriorPhi){
  theta.prime <- log(c(beta, phi))
  part1 <- sum(theta.prime*SUFF.BETA.PHI)
  JJ <- nrow(SUFF.MCMC)
  temp <- numeric(JJ)
  for (j in 1:JJ){
    temp[j] <- sum((theta.prime-THETA.FIXO)*SUFF.MCMC[j, ]) # < theta - psi, T*(N) >
  }
  part2 <- mean(exp(temp))
  like.theta <- (part1 - log(part2)) - gamma*SFALL[1] + log(gamma)*SFALL[3] + gamma*SFALL[4] -
    SFALL[1]*log(1 - exp(-gamma))
  return(logpriorPhi(phi, PHI.P) + like.theta)
}

posterioriGamma <- function(proc.fmat, beta, phi, gamma, sigma, kappa,
                            GAMMA.P = gamma.p,
                            # A1 = a1, A2 = a2, B1 = b1, B2 = b2,
                            # NN = N, Rc = R_centers, R = R_clusters,
                            SFALL = suff.all,
                            logpriorGamma){
  like.gamma <- -gamma*SFALL[1] + log(gamma)*SFALL[3] + gamma*SFALL[4] -
    SFALL[1]*log(1 - exp(-gamma))
  return(logpriorGamma(gamma, GAMMA.P) + like.gamma)
}

posterioriSigma <- function(proc.fmat, beta, phi, gamma, sigma, kappa,
                            SIGMA.P = sigma.p,
                            A1 = a1, A2 = a2, B1 = b1, B2 = b2,
                            NN = N, Rc = R_centers, R = R_clusters,
                            logpriorSigma){
  # return(log(dinvgamma(sigma, shape = sigma.p^2/10, scale = sigma.p)) +
  #!# MCMCpack; alpha (shape), beta (scale)
  suff.all <- sufficientStatAll(proc.fmat, sigma, kappa,
                                A1, A2, B1, B2, NN, Rc, R)
  like.mat <- -gamma*suff.all[1] + log(gamma)*suff.all[3] + gamma*suff.all[4] +
    suff.all[5] - suff.all[1]*log(1 - exp(-gamma))
  return(logpriorSigma(sigma, SIGMA.P) + like.mat)
}

posterioriKappa <- function(proc.fmat, beta, phi, gamma, sigma, kappa,
                            KAPPA.P = kappa.p,
                            A1 = a1, A2 = a2, B1 = b1, B2 = b2,
                            NN = N, Rc = R_centers, R = R_clusters,
                            logpriorKappa){
  suff.all <- sufficientStatAll(proc.fmat, sigma, kappa,
                                A1, A2, B1, B2, NN, Rc, R)
  like.mat <- -gamma*suff.all[1] + log(gamma)*suff.all[3] + gamma*suff.all[4] +
    suff.all[5] - suff.all[1]*log(1 - exp(-gamma))
  return(logpriorKappa(kappa, KAPPA.P) + like.mat)
}
