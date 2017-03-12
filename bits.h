#ifndef BMM_BITS_H
#define BMM_BITS_H

#include <stdbool.h>

inline unsigned char bmm_setb(unsigned char const x, unsigned char const n) {
  return (unsigned char) (x | (1 << n));
}

inline unsigned char bmm_clearb(unsigned char const x, unsigned char const n) {
  return (unsigned char) (x & ~(1 << n));
}

inline unsigned char bmm_complb(unsigned char const x, unsigned char const n) {
  return (unsigned char) (x ^ (1 << n));
}

inline void bmm_setbp(unsigned char* const x, unsigned char const n) {
  *x = bmm_setb(*x, n);
}

inline void bmm_clearbp(unsigned char* const x, unsigned char const n) {
  *x = bmm_clearb(*x, n);
}

inline void bmm_complbp(unsigned char* const x, unsigned char const n) {
  *x = bmm_complb(*x, n);
}

inline bool bmm_testb(unsigned char const x, unsigned char const n) {
  return ((x >> n) & 1) != 0;
}

#endif
