#include "sparsematrix.h"
#include "bitarray.h"


static int SparseMatrix_at_(const Dataset * const me, int n, int d);

void SparseMatrix_ctor(SparseMatrix* const me, int* p, int* i, int N, int D) {
  
  static struct DatasetVtbl const vtbl = {
    &SparseMatrix_at_
  };
  /* First call superclass' ctor */
  Dataset_ctor(&me->super, N, D);
  me->super.vptr = &vtbl;
  
  /* next, initialize the attributes specific to this subclass */ 
   
  Bitarray bitarray;
  Bitarray_ctor(&bitarray, p, i, N, D);
  me->_data = bitarray;
};

static int SparseMatrix_at_(Dataset const * const me, int n, int d) {
  
  SparseMatrix const * const m = (SparseMatrix const *) me; // Explicit downcast
  
  return Bitarray_at(&m->_data, n, d);
}
 