#ifndef BMM_ASET_H
/// Array-backed finite sets.
#define BMM_ASET_H

#include <stdbool.h>
#include <stddef.h>
#include <string.h>

#include "ext.h"
#include "size.h"

__attribute__ ((__nonnull__, __pure__))
inline bool bmm_aset_canins(size_t* const pnmemb, size_t const ncap) {
  return *pnmemb < ncap;
}

__attribute__ ((__nonnull__))
inline void bmm_aset_reins(void* restrict const paset,
    size_t* restrict const pnmemb, size_t const size,
    void const* restrict const memb) {
  unsigned char* const buf = paset;
  (void) memcpy(&buf[*pnmemb * size], memb, size);
}

__attribute__ ((__nonnull__))
inline void bmm_aset_ins(void* restrict const paset,
    size_t* restrict const pnmemb, size_t const size,
    void const* restrict const memb) {
  bmm_aset_reins(paset, pnmemb, size, memb);

  ++*pnmemb;
}

__attribute__ ((__nonnull__, __pure__))
inline bool bmm_aset_candel(size_t* const pnmemb, size_t const imemb) {
  return imemb < *pnmemb;
}

__attribute__ ((__nonnull__))
inline void bmm_aset_redel(void* restrict const paset,
    size_t* restrict const pnmemb, size_t const size,
    size_t const imemb) {
  unsigned char* const buf = paset;
  (void) memcpy(&buf[imemb * size], &buf[*pnmemb * size], size);
}

__attribute__ ((__nonnull__))
inline void bmm_aset_del(void* restrict const paset,
    size_t* restrict const pnmemb, size_t const size,
    size_t const imemb) {
  bmm_aset_redel(paset, pnmemb, size, imemb);

  --*pnmemb;
}

#define BMM_ASET_CANINS(aset, nmemb) \
  bmm_aset_canins(&nmemb, nmembof(aset))

#define BMM_ASET_REINS(aset, nmemb, memb) \
  bmm_aset_reins(aset, &nmemb, msizeof(aset), &memb)

#define BMM_ASET_INS(aset, nmemb, memb) \
  bmm_aset_ins(aset, &nmemb, msizeof(aset), &memb)

#define BMM_ASET_CANDEL(aset, nmemb) \
  bmm_aset_candel(&nmemb, nmembof(aset))

#define BMM_ASET_REDEL(aset, nmemb, imemb) \
  bmm_aset_redel(aset, &nmemb, msizeof(aset), imemb)

#define BMM_ASET_DEL(aset, nmemb, imemb) \
  bmm_aset_del(aset, &nmemb, msizeof(aset), imemb)

#endif
