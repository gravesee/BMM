#ifndef DATASET_H_
#define DATASET_H_

#include "stdlib.h"

struct DatasetVtbl; /* Forward declaration */

typedef struct _dataset {
  struct DatasetVtbl const *vptr; /* <== Datasets virtual pointer */
  
  // Base variables
  int N;
  int D; // number of rows and dimension

} Dataset;

struct DatasetVtbl {
  int (*at)(Dataset const * const me, int n, int d);
};

/* Datasets operations */

void Dataset_ctor(Dataset * const me, int N, int D);

// Generic operations on dataset inheritors
static inline int at(Dataset * const me, int n, int d) {
  return (*me->vptr->at)(me, n, d);
}

#endif /* DATASET_H_ */

