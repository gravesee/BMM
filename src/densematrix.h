#ifndef DENSEMATRIX_H
#define DENSEMATRIX_H

#include "assert.h"
#include "R.h"
#include "Rinternals.h"
#include "dataset.h"

typedef struct _densematrix {
  Dataset super;
  int* _data;
} DenseMatrix;

void DenseMatrix_ctor(DenseMatrix* const me, SEXP matrix, int N, int D);

#endif /* DENSEMATRIX_H */
