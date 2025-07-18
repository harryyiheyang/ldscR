#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::rowvec vectorMatrixDotProduct(const arma::vec& v, const arma::mat& M) {
  if (v.size() != M.n_rows) {
    stop("The vector length and the number of rows in the matrix must be the same.");
  }

  arma::rowvec results = v.t() * M;

  return results;
}
