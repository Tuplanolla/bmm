// Bitwise operations.
#ifndef BMM_BIT_H
#define BMM_BIT_H

#include <stdbool.h>

#include "ext.h"

// The call `bmm_bit_test(x, n)`
// checks whether `x` has the bit at index `n` set.
__attribute__ ((__const__, __pure__))
inline bool bmm_bit_test(unsigned char const x, unsigned char const n) {
  return ((x >> n) & 1) != 0;
}

// The call `bmm_bit_set(x, n)`
// returns `x` with the bit at index `n` set.
__attribute__ ((__const__, __pure__))
inline unsigned char bmm_bit_set(unsigned char const x,
    unsigned char const n) {
  return (unsigned char) (x | (1 << n));
}

// The call `bmm_bit_clear(x, n)`
// returns `x` with the bit at index `n` unset.
__attribute__ ((__const__, __pure__))
inline unsigned char bmm_bit_clear(unsigned char const x,
    unsigned char const n) {
  return (unsigned char) (x & ~(1 << n));
}

// The call `bmm_bit_compl(x, n)`
// returns `x` with the bit at index `n` flipped.
__attribute__ ((__const__, __pure__))
inline unsigned char bmm_bit_compl(unsigned char const x,
    unsigned char const n) {
  return (unsigned char) (x ^ (1 << n));
}

// The call `bmm_bit_pset(ptr, n)`
// sets the bit at index `n` in `ptr`.
__attribute__ ((__nonnull__))
inline void bmm_bit_pset(unsigned char* const ptr, unsigned char const n) {
  *ptr = bmm_bit_set(*ptr, n);
}

// The call `bmm_bit_pclear(ptr, n)`
// unsets the bit at index `n` in `ptr`.
__attribute__ ((__nonnull__))
inline void bmm_bit_pclear(unsigned char* const ptr, unsigned char const n) {
  *ptr = bmm_bit_clear(*ptr, n);
}

// The call `bmm_bit_pcompl(ptr, n)`
// flips the bit at index `n` in `ptr`.
__attribute__ ((__nonnull__))
inline void bmm_bit_pcompl(unsigned char* const ptr, unsigned char const n) {
  *ptr = bmm_bit_compl(*ptr, n);
}

#endif
