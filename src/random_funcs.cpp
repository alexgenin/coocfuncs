
#include <RcppArmadillo.h>

using namespace Rcpp; 

// Generate random integer between min and max
int randn(double min, double max) { 
  NumericVector result = runif(1, 0, 1);
  return floor( min + result[0] * ((max+1) - min) );
}

// Generate random number between 0 and 1 as scalar using runif function in 
// Rcpp. Chapuza.
// [[Rcpp::export]] // for debug purposes only
double randp() { 
  NumericVector result = runif(1,0,1);
  return result[0];
}
 
