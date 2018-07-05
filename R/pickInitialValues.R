#' @export
pickInitialValues <- function(){
  list(iniBeta  = function(fmat, win) length(fmat)/area.owin(win),
       iniPhi   = function() 2,
       iniGamma = function(fmat) max(sapply(fmat, function(z) length(z$x))),
       iniSigma = function(fmat) sqrt(sum(sapply(fmat, function(z) {
         if(length(z$x) > 0){
           sum((z$x-z$centers[1])^2 + (z$y-z$centers[2])^2)
         }} )) / sum(sapply(fmat, function(z) length(z$x)))),
       iniKappa = function(angles) mle.vonmises(circular(angles, modulo = "2pi"), mu = circular(0))$kappa)
}
