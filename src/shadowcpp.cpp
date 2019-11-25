#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double shadowcpp(const int N, const arma::mat& X,
                 const int Cc,
                 const arma::vec& x,
                 const arma::vec& y,
                 const arma::vec& t,
                 const double sigma, const double r) {

  double area_shadow;
  double temp = 0;
  arma::mat d(N, Cc);
  arma::vec a(Cc);
  arma::vec temp1(N);
  d.zeros(N, Cc);
  a.zeros(Cc);
  temp1.zeros(N);

  for(int ic = 0; ic < Cc; ic++) { // For each finger
    for(int jc = 0; jc < N; jc++) { // For each Monte Carlo shadow evaluation...
      //      printf("jc %i \n", jc);
      temp = (X(jc,0)-x(ic))*(X(jc,0)-x(ic)) + (X(jc,1)-y(ic))*(X(jc,1)-y(ic));
      if(temp <= r*r){
        d(jc,ic) = 1; // Candidate is within shadow of given center
      } else {
        d(jc,ic) = 0; // Candidate is outside shadow of given center
      }
    }
  }

  for(int ij = 0; ij < N; ij++) {
    temp1(ij) = d(ij,0); // Within R distance from the first finger
  }

  double atemp = mean(temp1); // proportion near first finger
  a(0) = atemp*(1-t(0)); // ConeArea > t, since t \in [0,1]
  //        printf("0, a(0) %f \n", a(0,0));

  //!/ if (Cc > 1){ // TRUE by construction

  for(int ic = 1; ic < Cc; ic++){ // For each finger, starting at finger 2...
    temp1.ones(N); // Reset temp1
    for(int ii = 0; ii < ic; ii++){ // For each previous finger
      for(int ij = 0; ij < N; ij++){
        temp1(ij) *= (1-d(ij,ii)); // if it was inside any previous finger area, set to zero
      }
    }
    for(int ij = 0; ij < N; ij++){
      temp1(ij) *= d(ij,ic); // if it is inside current finger area, maintain at 1
    }
    atemp = mean(temp1); // proportion inside current finger area, exclusively
    a(ic) = atemp*(1-t(ic)); // current finger ConeArea > t, since t \in [0,1]
    //      printf("%i, a(ic) %f \n", ic, a(ic,0));
  }
  //!/ }

  area_shadow = sum(a);

  return area_shadow;
}
