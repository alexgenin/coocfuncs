// 
// 
// This function computes the absolute value of overlap between two segments
// 

#include <RcppArmadillo.h>

using namespace Rcpp; 


// Define here the association function you want to use
#define ASSOCF overlap_binary
// #define ASSOCF overlap_abs

// Define the tolerance (number in cm by which plants are grown) when 
//   considering overlapping
#define TOL 1

#define MIN(a,b) a<b ? a : b

//
// 
// Association functions !
// 

// Absolute length of overlap of species i and j
//[[Rcpp::export]]
double overlap(double xi_i, 
               double xi_j, 
               double xe_i, 
               double xe_j) { 
  
  // Case where the individual is the same (account for floating-point errors)
  if ( std::fabs(xi_j - xi_i) < .00001 && std::fabs(xe_j - xe_i) < .00001 ) { 
    return( xe_i - xi_i );
  }
//   Rcout << "i: " << xi_i << "-" << xe_i << " and j: " << xi_j << "-" << xe_j << ": "; 
  
  double ans = 0; 
  
  // Case where there is no overlap 
  if ( xe_i < xi_j || xe_j < xi_i ) { 
    ans = 0; 
  }
  
  // Case where i is included in j -> we return the length of i
  if ( xi_j <= xi_i && xe_i <= xe_j ) { 
    ans = MIN(xe_i - xi_i, xe_j - xi_j);
  }
  
  // Case where j is included in i -> we return the length of j
  if ( xi_i <= xi_j && xe_j <= xe_i ) {
    ans = MIN(xe_i - xi_i, xe_j - xi_j);
  }
  
  // Case where i is slightly after j 
  if ( xi_i <= xe_j && xi_j <= xi_i ) {
    ans = xe_j - xi_i; 
  }
  
  //  Case where i is slightly before j 
  if ( xi_j <= xe_i && xe_i <= xe_j) { 
    ans = xe_i - xi_j;
  }
  
//   Rcout << ans <<  "\n"; 
  
  return(ans); // no overlap
}


// Binary value whether there is overlap or not
//[[Rcpp::export]]
double overlap_binary(double istart, 
                      double iend, 
                      double jstart, 
                      double jend) { 
  
  // Case  -----(+++i++)-----|+j++|------
  if ( (istart-TOL) > (jend+TOL) || (iend+TOL) < (jstart-TOL) ) { 
    return(0.0);
  } else { 
  // Case  -----(+++i++)-----------------
  //       --------|+j++|----------------
    return(1.0);
  }
  
}


