#include "dataset.h"
#include "R.h"


typedef struct _sparsematrix {
  Dataset base;
  int* _p; // column pointer
  int* _i; // row index
} SparseMatrix;

int SparseMatrixAt(void* sparseMatrix, int n, int d) {
  
  SparseMatrix* s = (SparseMatrix*) sparseMatrix;
  
  // bounds check
  //assert(n <= m->base.N && n >= 0);
  //assert(d <= m->base.D && d >= 0);
  
  Rprintf("Must implement\n");
  
  return 0;
  
  //return m->_data[d * m->base.N + n];
  
}

SparseMatrix* getSparseMatrixInstance(SEXP p, SEXP i, int N, int D) {
  SparseMatrix* m = (SparseMatrix*)malloc(sizeof(SparseMatrix));
  m->_p = INTEGER(p);
  m->_i = INTEGER(i);
  m->base.N = N;
  m->base.D = D;
  m->base.functions.at = SparseMatrixAt;
  return m;
}
