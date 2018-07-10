#include "assert.h"
#include "R.h"
#include "Rinternals.h"
#include "densematrix.h"

static int DenseMatrix_at_(const Dataset * const me, int n, int d);

/* constructor implementation */
void DenseMatrix_ctor(DenseMatrix* const me, SEXP matrix, int N, int D) {
  static struct DatasetVtbl const vtbl = {
    &DenseMatrix_at_
  };
  /* First call superclass' ctor */
  Dataset_ctor(&me->super, N, D);
  me->super.vptr = &vtbl;
  
  /* next, initialize the attributes specific to this subclass */
  me->_data = INTEGER(matrix);
  
}

static int DenseMatrix_at_(Dataset const * const me, int n, int d) {
  
  DenseMatrix const * const m = (DenseMatrix const *) me; // Explicit downcast
  
  // bounds check
  assert(n <= m->super.N && n >= 0);
  assert(d <= m->super.D && d >= 0);
  
  return m->_data[d * m->super.N + n];
  
}
