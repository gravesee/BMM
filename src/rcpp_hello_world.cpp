
#include <Rcpp.h>
#include "densematrix.h"
using namespace Rcpp;

// [[Rcpp::export]]
SEXP test_dense_matrix(SEXP m, int n, int d) {

  Dataset* matrix = (Dataset*) getDenseMatrixInstance(m, n, d);
  
  Rprintf("(1,1) = (%d)\n", matrix->functions.at(matrix, 0, 0));
  Rprintf("(2,2) = (%d)\n", matrix->functions.at(matrix, 1, 1));
  Rprintf("(1,2) = (%d)\n", matrix->functions.at(matrix, 0, 1));
  Rprintf("(2,1) = (%d)\n", matrix->functions.at(matrix, 1, 0));
  
  // Method for freeing it
  free(matrix);
  
  return R_NilValue;
}
