#ifndef BMM_ASTACK_H
/// Array-backed finite stacks.
#define BMM_ASTACK_H

#include <stdbool.h>
#include <stddef.h>
#include <string.h>

#include "ext.h"
#include "size.h"

__attribute__ ((__nonnull__, __pure__))
inline bool bmm_astack_canins(size_t* const pnmemb, size_t const ncap) {
  return *pnmemb < ncap;
}

__attribute__ ((__nonnull__))
inline void bmm_astack_preins(void* restrict const pastack,
    size_t* restrict const pnmemb, size_t const size,
    void const* restrict const memb) {
  unsigned char* const buf = pastack;
  (void) memcpy(&buf[*pnmemb * size], memb, size);
}

__attribute__ ((__nonnull__))
inline void bmm_astack_ins(void* restrict const pastack,
    size_t* restrict const pnmemb, size_t const size,
    void const* restrict const memb) {
  bmm_astack_preins(pastack, pnmemb, size, memb);

  ++*pnmemb;
}

__attribute__ ((__nonnull__, __pure__))
inline bool bmm_astack_candel(size_t* const pnmemb) {
  return *pnmemb > 0;
}

__attribute__ ((__nonnull__))
inline void bmm_astack_del(size_t* const pnmemb) {
  --*pnmemb;
}

#define BMM_ASTACK_CANINS(astack, nmemb) \
  bmm_astack_canins(&nmemb, nmembof(astack))

#define BMM_ASTACK_PREINS(astack, nmemb, memb) \
  bmm_astack_preins(astack, &nmemb, msizeof(astack), &memb)

#define BMM_ASTACK_INS(astack, nmemb, memb) \
  bmm_astack_ins(astack, &nmemb, msizeof(astack), &memb)

#define BMM_ASTACK_CANDEL(astack, nmemb) \
  bmm_astack_candel(&nmemb, nmembof(astack))

#define BMM_ASTACK_DEL(astack, nmemb) \
  bmm_astack_del(&nmemb)

#endif
