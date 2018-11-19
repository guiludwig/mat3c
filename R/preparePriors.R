#' @export
preparePriors <- function(){
  list(priorBeta = function(beta, beta.p) dlnorm(beta, log(beta.p), 1, log = TRUE),
       priorPhi = function(phi, phi.p) {
         if(phi > 1){
           # Divide area by 1-alpha so it integrates to 1
           log(0.8) + dgamma(phi, shape = 18/3, scale = 1/3, log = TRUE) - pgamma(1, shape = 18/3, scale = 1/3, lower.tail = FALSE, log = TRUE)
         } else if(phi == 1) {
           log(0.2)
         } else {
           -Inf
         }
       },
       priorGamma = function(gamma, gamma.p) dgamma(gamma, shape = gamma.p, rate = 1, log = TRUE),
       priorSigma = function(sigma, sigma.p) ((sigma.p^2)/10) * log(sigma.p) - lgamma((sigma.p^2)/10) - ((sigma.p^2)/10 + 1) * log(sigma) - (sigma.p/sigma),
       priorKappa = function(kappa, kappa.p) dgamma(kappa, shape = (kappa.p^2)/10, scale = 10/kappa.p, log = TRUE))
}
