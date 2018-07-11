#ifndef BITARRAY_H
#define BITARRAY_H
#include <limits.h>

#define BITS_IN_INT (sizeof(int) * 8)

/* Hat-tip to SO for these macros:
* "How to define and work with an array of bits in C?"
* https://stackoverflow.com/a/30590727/919872
*/

#define SetBit(A,k)     ( A[(k / BITS_IN_INT)] |= (1 << (k % BITS_IN_INT)) )
#define TestBit(A,k)    ( A[(k / BITS_IN_INT)] & (1 << (k % BITS_IN_INT)) )

/*
 * Bitarray structure for storing binary data as array of ints
 * \param int** data an array of array of ints
 * \param nrow int storing number of elements to store -- usually a "row" of data from a sparse matrix
 * \param ncol int storing number of bytes needed to store all of the bits
 * \param nbits int storing number of individual bits to represent
 */
typedef struct  {
  int** _data;
  int _nrows;
  int _nbits;
} Bitarray;


int Bitarray_at(const Bitarray * me, int n, int d);

void Bitarray_set(Bitarray * me, int n, int d);

/*
 * Helper method returning number of bytes need to store nbits for each machine
 * Calculated as the # of bits needed / # bits in an int
 */
int bytes_needed(int d);

/*
 * Allocate and zero out the space needed to store requested elements
 */
int** allocate_bytes(int n, int d);

/*
 * Constructor initializing a Bitarray struct and allocating spaces
 */
void Bitarray_ctor(Bitarray * const me, int*p, int*i, int n, int d);


/*
 * Destructor freeing memory allocated for Bitarray
 */
void Bitarray_free(Bitarray * const me);


#endif /* BITARRAY_H */


