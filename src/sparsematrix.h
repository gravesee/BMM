#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include "R.h"
#include "Rinternals.h"
#include "dataset.h"

typedef struct _sparsematrix {
  Dataset super;
  int* _p;
  int* _i;
} SparseMatrix;

void SparseMatrix_ctor(SparseMatrix* const me, SEXP p, SEXP i);

#endif /* SPARSEMATRIX_H */
