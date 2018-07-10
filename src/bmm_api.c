
#include "R.h"
#include "Rinternals.h"
#include "densematrix.h"
#include "bmm_em.h"

SEXP test_dense_matrix(SEXP m, SEXP n, SEXP d) {
   
  DenseMatrix x;
  DenseMatrix_ctor(&x, m, asInteger(n), asInteger(d));
  
  Rprintf("(1,1) = (%d)\n", at(&x.super, 0,0));
  Rprintf("(1,2) = (%d)\n", at(&x.super, 1,0));
  Rprintf("(2,1) = (%d)\n", at(&x.super, 0,1));
  Rprintf("(2,2) = (%d)\n", at(&x.super, 1,1));
  
  
  return R_NilValue;
}

SEXP bmm_dense_matrix(SEXP m, SEXP n, SEXP d, SEXP K, SEXP max_iter, SEXP verbose) {
  
  DenseMatrix x;
  DenseMatrix_ctor(&x, m, asInteger(n), asInteger(d));
  
  bmm_em_result res = em(&x.super, asInteger(K), asInteger(max_iter), asInteger(verbose));
  
  // Pass pointer to keep track of R protect uses
  // Convert to SEXP to return to R
  int prtCnt = 0;
  SEXP out = convert_bmm_em_result(&x.super, &res, &prtCnt);
  
  free_bmm_em_result(&res);
  
  UNPROTECT(prtCnt);
  
  return out;
} 