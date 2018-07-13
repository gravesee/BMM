
#include "R.h"
#include "Rinternals.h"
#include "densematrix.h"
#include "sparsematrix.h"
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



SEXP test_sparse_matrix(SEXP p, SEXP i, SEXP n, SEXP d) {
  
  SparseMatrix x;
  SparseMatrix_ctor(&x, INTEGER(p), INTEGER(i), asInteger(n), asInteger(d));
  
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



SEXP bmm_sparse_matrix(SEXP p, SEXP i, SEXP n, SEXP d, SEXP K, SEXP max_iter, SEXP verbose) {
  
  SparseMatrix x;
  SparseMatrix_ctor(&x, INTEGER(p), INTEGER(i), asInteger(n), asInteger(d));
  
  bmm_em_result res = em(&x.super, asInteger(K), asInteger(max_iter), asInteger(verbose));
  
  Bitarray_free(&x._data); // free data holding bitarray
  
  // Pass pointer to keep track of R protect uses
  // Convert to SEXP to return to R
  int prtCnt = 0;
  SEXP out = convert_bmm_em_result(&x.super, &res, &prtCnt);
  
  free_bmm_em_result(&res);
  
  UNPROTECT(prtCnt);
  
  return out;
} 


SEXP predict_dense_matrix(SEXP m, SEXP protos, SEXP pis, SEXP n, SEXP d, SEXP K) {
  
  int _K = asInteger(K);
  int _d = asInteger(d);
  int _n = asInteger(n);
  
  DenseMatrix x;
  DenseMatrix_ctor(&x, m, _n, _d);
  
  // Must pass transposed protos from R!
  double ** _protos = (double **) malloc(_K * _d * sizeof(double*));
  for (int k = 0; k < _K; k++) {
    _protos[k] = &REAL(protos)[k * _d];
  }
  
  znk_result res = predict_log_z_nk(&x.super, REAL(pis), _protos, _K);
  
  // Pass pointer to keep track of R protect uses
  // Convert to SEXP to return to R
  int prtCnt = 0;
  SEXP out = convert_znk_result(&x.super, &res, &prtCnt);
  
  free_znk_result(&res);
  free(_protos);
  
  UNPROTECT(prtCnt);
  
  return out;
} 


SEXP predict_sparse_matrix(SEXP p, SEXP i, SEXP protos, SEXP pis, SEXP n, SEXP d, SEXP K) {
  
  int _K = asInteger(K);
  int _n = asInteger(n);
  int _d = asInteger(d);
  
  SparseMatrix x;
  SparseMatrix_ctor(&x, INTEGER(p), INTEGER(i), _n, _d);
  
  // Must pass transposed protos from R!
  double ** _protos = (double **) malloc(_K * _d * sizeof(double*));
  for (int k = 0; k < _K; k++) {
    _protos[k] = &REAL(protos)[k * _d];
  }
  
  znk_result res = predict_log_z_nk(&x.super, REAL(pis), _protos, _K);
  
  // Pass pointer to keep track of R protect uses
  // Convert to SEXP to return to R
  int prtCnt = 0;
  SEXP out = convert_znk_result(&x.super, &res, &prtCnt);
  
  free_znk_result(&res);
  free(_protos);
  
  UNPROTECT(prtCnt);
  
  return out;
} 



