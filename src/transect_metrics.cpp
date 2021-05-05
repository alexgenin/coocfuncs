// 
// Compute the vegetation cover in a transect
// 

#include <RcppArmadillo.h>
using namespace arma; 
using namespace Rcpp; 

// Defined in shuffle_core.cpp
double get_cover(arma::mat& xs, double xmax); 

//[[Rcpp::export]]
double trans_cover(arma::vec& xis, 
                   arma::vec& xes, 
                   double xmax) { 
  
  arma::mat xs(xis.n_elem, 2); 
  xs.col(0) = xis; 
  xs.col(1) = xes; 
  
  double cover = get_cover(xs, xmax); 
  
  return( cover ); 
}

