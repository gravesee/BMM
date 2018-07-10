#include "dataset.h"
#include "R.h"

typedef struct _densematrix {
  Dataset base;
  int* _data;
} DenseMatrix;

int DenseMatrixAt(void* denseMatrix, int n, int d) {
  
  DenseMatrix* m = (DenseMatrix*) denseMatrix;
  // bounds check
  assert(n <= m->base.N && n >= 0);
  assert(d <= m->base.D && d >= 0);
  
  return m->_data[d * m->base.N + n];
  
}

DenseMatrix* getDenseMatrixInstance(SEXP matrix, int N, int D) {
  DenseMatrix* m = (DenseMatrix*)malloc(sizeof(DenseMatrix));
  m->_data = INTEGER(matrix);
  m->base.N = N;
  m->base.D = D;
  m->base.functions.at = DenseMatrixAt;
  return m;
}
