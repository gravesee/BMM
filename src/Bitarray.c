#include "bitarray.h"
#include "stdlib.h"

void Bitarray_ctor(Bitarray * const me, int*p, int*i, int n, int d) {
  me->_data = allocate_bytes(n, d);
  me->_nrows = n;
  me->_nbits = d;
  
  // Fill _data with integers that correspond to set bits
  
  // loop over pointer indices
  for (int bit = 0; bit < d; bit++) {
    int start = p[bit];
    int stop = p[bit + 1];
    // loop over rows
    for (int r = start; r < stop; r++) {
      int row = i[r];
      Bitarray_set(me, row, bit);
    }
  }
  
};

int bytes_needed(int nbits) {
  return nbits / BITS_IN_INT + 1;
}

void Bitarray_free(Bitarray * const me) {
  for (int i = 0; i < me->_nrows; i++) {
    free(me->_data[i]);
  }
  free(me->_data);
}

int** allocate_bytes(int n, int d) {
  int ncol = bytes_needed(d);
  
  int** m = (int**) calloc(n, sizeof(int*));
  
  // allocate arrays of ints representing each row
  for (int r = 0; r < n; r++) {
    m[r] = (int*) calloc(ncol, sizeof(int));
  }
  
  return m;
}


int Bitarray_at(const Bitarray * me, int n, int d) {
  return (int) (TestBit(me->_data[n], d) != 0);
};


void Bitarray_set(Bitarray * const me, int n, int d) {
  SetBit(me->_data[n], d);
};

