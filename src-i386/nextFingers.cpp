#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::LogicalVector nextFingers(double markedN, Rcpp::NumericMatrix fingers0,
                                double Rclusters) {
  Rcpp::LogicalVector indFingers(markedN);
  for(int i = 0; i < markedN; i++) indFingers[i] = true;
  bool test;
  if(markedN > 1){ // At least two fingers...
    int tempN = fingers0.nrow();
    for(int i = 0; i < tempN-1; i++){
      if(indFingers[i]){
        for(int j = i+1; j < tempN; j++){
          test = ((fingers0(i,0) - fingers0(j,0))*(fingers0(i,0) - fingers0(j,0)) + (fingers0(i,1) - fingers0(j,1))*(fingers0(i,1) - fingers0(j,1)) > Rclusters*Rclusters);
          indFingers[j] = indFingers[j] && test;
        }
      }
    }
  }
  return(indFingers);
}

/*** R
set.seed(1)
n <- 100 # n <- 50000
fingers0 <- matrix(rnorm(2*n), ncol=2)
RC <- .5
system.time({
  a <- nextFingers(n, fingers0, RC)
})

system.time({
  ind.fingers <- rep(TRUE, n)
  #!# Only allow fingers that are far from existing ones to be born
  #!# Note: A lot of fingers are born (~50), only 4-5 survive in general
  if (n > 1){ # At least two fingers...
    tempN <- nrow(fingers0)
    for(i in 1:(tempN-1)){
      if(ind.fingers[i]){ #!# SLOW
        ind.fingers[(i + 1):tempN] <- ind.fingers[(i + 1):tempN] & ((fingers0[i, 1] - fingers0[(i + 1):tempN, 1])^2 +
                                                                      (fingers0[i, 2] - fingers0[(i + 1):tempN, 2])^2 > RC^2)
      }
    }
  }
})
all.equal(a, ind.fingers)
*/
