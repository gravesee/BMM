#include "dataset.h"
#include "assert.h"

static int Dataset_at_(Dataset const * const me, int n, int d);

void Dataset_ctor(Dataset * const me, int N, int D) {
  static struct DatasetVtbl const vtbl = {
    &Dataset_at_
  };
  
  me->vptr = &vtbl; /* Hook the vptr to the vtbl */
  me->N = N;
  me->D = D;
};

static int Dataset_at_(Dataset const * const me, int n, int d) {
  assert(0);
  return 0;
};

