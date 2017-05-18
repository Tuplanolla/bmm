#ifndef BMM_ASET_H
/// Array-backed finite sets.
#define BMM_ASET_H

#include <stdbool.h>
#include <stddef.h>
#include <string.h>

#include "ext.h"
#include "size.h"

inline bool bmm_aset_canins(size_t* const pnmemb, size_t const ncap) {
  return *pnmemb < ncap;
}

inline void bmm_aset_ins(void* restrict const paset,
    size_t* restrict const pnmemb, size_t const size,
    void const* restrict const memb) {
  size_t const nmemb = *pnmemb;

  unsigned char* const buf = paset;
  (void) memcpy(&buf[nmemb * size], memb, size);

  *pnmemb = nmemb + 1;
}

inline bool bmm_aset_candel(size_t* const pnmemb, size_t const imemb) {
  return imemb < *pnmemb;
}

inline void bmm_aset_del(void* restrict const paset,
    size_t* restrict const pnmemb, size_t const size,
    size_t const imemb) {
  size_t const nmemb = *pnmemb;

  unsigned char* const buf = paset;
  (void) memcpy(&buf[imemb * size], &buf[nmemb * size], size);

  *pnmemb = nmemb - 1;
}

#define BMM_ASET_CANINS(aset, nmemb) \
  bmm_aset_canins(&nmemb, nmembof(aset))

#define BMM_ASET_INS(aset, nmemb, memb) \
  bmm_aset_ins(aset, &nmemb, msizeof(aset), memb)

#define BMM_ASET_CANDEL(aset, nmemb) \
  bmm_aset_candel(&nmemb, nmembof(aset))

#define BMM_ASET_DEL(aset, nmemb, imemb) \
  bmm_aset_del(aset, &nmemb, msizeof(aset), imemb)

#endif
