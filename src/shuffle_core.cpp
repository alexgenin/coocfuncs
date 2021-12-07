// 
// This function takes a transect and shuffles its elements
// 

#include <RcppArmadillo.h>
#include "shared_functions.h"
#include <stdlib.h>

using namespace Rcpp; 
using namespace arma; 


arma::mat shuffle_inplace(arma::mat xs, 
                          double xmax) { 
  
  uword Ninds = xs.n_rows; 
  
  for (uword i=0; i<Ninds; i++) { 
    // Build new coordinates
    double indsize = xs(i, 1) - xs(i, 0); 
    double newpos = fmod(xmax + as_scalar(randu(1)) * (xmax - indsize), xmax); 
    
    xs(i, 0) = newpos; 
    xs(i, 1) = newpos + indsize; 
  }
  
  return(xs); 
}

//[[Rcpp::export]]
arma::mat shuffle_core(arma::mat& xs, double xmax) { 
  
  arma::mat newxs = shuffle_inplace(xs, xmax); 
  
  return(newxs); 
}

//'@export
//[[Rcpp::export]]
double get_cover(arma::mat& xs, 
                 double xmax) { 

  // Step along the transect. We choose either to take a set number of points, 
  // or we take the minimum size divided by x
  double xstep = min(xs.col(1) - xs.col(0)) / 1.5; 
  double xstep2 = xmax / (1 + 1000); 
  xstep = xstep < xstep2 ? xstep : xstep2;
  
  // Number of positions with vegetation along the transect
  int ntotal = 0; 
  int ntested = 0; 
  
  // We are going to go along indivs in sorted order. This way we know 
  // we can only check indivs after a current position and check less things 
  // this way. 
  arma::uvec ind_order = sort_index(xs.col(1)); 
  uword minind_sorted = 0; 
  
  // Go along the transect
  // x is current position along the transect 
  for (double x=0; x<=xmax; x+=xstep) { 
    
    bool someone_in_x = false; 
    uword ind_sort = minind_sorted; 
    while ( ! someone_in_x && ind_sort < ind_order.n_elem ) { 
      // If the position is after the end of the indiv, then discard this 
      // indiv from next iter checks. 
//       Rcout << "indsort: " << ind_sort << " -> " << ind_order(ind_sort) << "\n"; 
      uword ind = ind_order(ind_sort); 
      if ( x > xs(ind, 1) ) { 
        minind_sorted = ind_sort; 
//         Rcout << "ind " << ind << "is below x: " << x << "\n"; 
      // Else we test if the position is occupied i.e. x is after the beginning 
      // of the individual. 
      } else if ( x >= xs(ind, 0) ) { 
        someone_in_x = true; 
        ntotal++; 
      }
      ind_sort++; 
    }
    
    ntested++; 
//     Rcout << x << " by " << xstep << ":" << ntotal << "/" << ntested << "\n"; 
  }
  
  return( (double)(ntotal) / ntested ); 
}
