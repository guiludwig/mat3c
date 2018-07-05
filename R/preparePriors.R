#' @export
preparePriors <- function(){
  list(priorBeta = function(beta, beta.p) dlnorm(beta, log(beta.p), 1, log = TRUE),
       priorPhi = function(phi, phi.p) dgamma(phi, shape = 8/3, scale = 3/4, log = TRUE),
       priorGamma = function(gamma, gamma.p) dgamma(gamma, shape = gamma.p, rate = 1, log = TRUE),
       priorSigma = function(sigma, sigma.p) ((sigma.p^2)/10) * log(sigma.p) - lgamma((sigma.p^2)/10) - ((sigma.p^2)/10 + 1) * log(sigma) - (sigma.p/sigma),
       priorKappa = function(kappa, kappa.p) dgamma(kappa, shape = (kappa.p^2)/10, scale = 10/kappa.p, log = TRUE))
}
