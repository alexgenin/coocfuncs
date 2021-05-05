// 
// 
// This file computes overlap between species in a transect
//

#include <RcppArmadillo.h>

#include "shared_functions.h"

using namespace Rcpp; 
using namespace arma; 

//[[Rcpp::export]]
arma::mat cont_overlap_core(arma::uword N, 
                            arma::uword Nsp, 
                            arma::vec& xis, 
                            arma::vec& xes, 
                            arma::uvec& attributs, 
                            arma::uvec& transectid, 
                            double tol) {  // as integers starting from 0
  
  arma::mat ovlps(Nsp, Nsp); 
  ovlps.fill(0); 
  
  for (uword i=0; i < N; i++) { 
    uword spi = attributs(i); 
    
    for (uword j=i+1; j < N; j++) {  // starts from i+1 to skip self-interaction
      uword spj = attributs(j); 
      
      double overlap_current = overlap(xis(i) - tol, xis(j), 
                                       xes(i) + tol, xes(j)); 
//    Rcout << "current overlap: " << overlap_current << "\n"; 
      
      // Check if the plants belong to the same transect
      if ( transectid(i) == transectid(j) ) { 
        ovlps(spi, spj) += overlap_current; 
        ovlps(spj, spi) += overlap_current; 
      }
    }
  }
  
  return(ovlps);
}
