#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include "R.h"
#include "Rinternals.h"
#include "dataset.h"
#include "bitarray.h"


typedef struct _sparsematrix {
  Dataset super;
  Bitarray _data;
} SparseMatrix;

void SparseMatrix_ctor(SparseMatrix* const me, int* p, int* i, int N, int D);

#endif /* SPARSEMATRIX_H */
