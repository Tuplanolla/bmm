// Bitwise operations.
#ifndef BMM_BIT_H
#define BMM_BIT_H

#include <stdbool.h>

inline unsigned char bmm_bit_set(unsigned char const x,
    unsigned char const n) {
  return (unsigned char) (x | (1 << n));
}

inline unsigned char bmm_bit_clear(unsigned char const x,
    unsigned char const n) {
  return (unsigned char) (x & ~(1 << n));
}

inline unsigned char bmm_bit_compl(unsigned char const x,
    unsigned char const n) {
  return (unsigned char) (x ^ (1 << n));
}

inline void bmm_bit_pset(unsigned char* const x, unsigned char const n) {
  *x = bmm_bit_set(*x, n);
}

inline void bmm_bit_pclear(unsigned char* const x, unsigned char const n) {
  *x = bmm_bit_clear(*x, n);
}

inline void bmm_bit_pcompl(unsigned char* const x, unsigned char const n) {
  *x = bmm_bit_compl(*x, n);
}

inline bool bmm_bit_test(unsigned char const x, unsigned char const n) {
  return ((x >> n) & 1) != 0;
}

#endif
