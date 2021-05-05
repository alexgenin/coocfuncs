//
//
// Helper function that bins data from a continuous df
//

#include "shared_functions.h"

#include <RcppArmadillo.h>

using namespace arma;

//[[Rcpp::export]]
arma::umat bin_data_core(arma::vec xis,
                         arma::vec xes,
                         arma::uvec spp,
                         arma::vec bins,
                         double binsize) {

  int n = bins.n_elem;
  int nsp = max(spp);
  int nindiv = xis.n_elem;

  arma::umat result(n, nsp);
  result.fill(0);

  // Loop over individuals
  for (int i=0; i<nindiv; i++) {
    // Loop over bins
    for (int b=0; b<n; b++) {
      double icenter   = (xes(i) + xis(i)) / 2;
      double ihalfsize = (xes(i) - xis(i)) / 2;

      if ( std::fabs(bins(b) - icenter) < ihalfsize ) {
        // Factors start at 1 in R so we need to adjust the index
        result(b, spp(i)-1) += 1;
      }
    }
  }

  return(result);
}

//[[Rcpp::export]]
arma::mat sample_transect_core(int N,
                                arma::vec xis,
                                arma::vec xes,
                                arma::uvec spp,
                                arma::uvec groups,
                                double xmax,
                                double binsize) {

  int nsp = max(spp);
  int ngroups = max(groups);
  uword nindiv = xis.n_elem;

  arma::mat result(N, nsp+2);
  result.fill(0);

  int s = 0;
  while (s < N) {

    // Choose a random bin
    double rbin = xmax;

    rbin = randp() * (xmax-binsize);

    // Choose a random group
    uword rgrp = randn(1, (double)(ngroups));

    // Check who is in there

    // Loop over individuals
    for (uword i=0; i<nindiv; i++) {

      // If in the right transect
      if ( groups(i) == rgrp ) {
//         double icenter   = (xes(i) + xis(i)) / 2;
//         double ihalfsize = (xes(i) - xis(i)) / 2;
//         double bincenter = (rbin + binsize) / 2;

        if ( overlap(rbin, xis(i), rbin+binsize, xes(i)) > 0 ) {
          // Factors start at 1 in R so we need to adjust the spp index
          result(s, 2+spp(i)-1) += 1.0;
        }

      }
    }

    // Save group and bin info
    result(s, 0) = rgrp;
    result(s, 1) = rbin;
    s++;

  }

  return(result);
}

