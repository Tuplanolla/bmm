#ifndef BMM_SORT_H
/// Higher-order sorting procedures.
#define BMM_SORT_H

#include <stddef.h>

#include "ext.h"

/// The call `bmm_sort_heapsort(nmemb, compar, swap, cls)`
/// uses unstable heap sort to rearrange `nmemb` items
/// with the comparison function `compar` and the swap procedure `swap`.
/// The closure `cls` is passed through to `compar` and `swap`.
/// Overflows are handled appropriately.
__attribute__ ((__nonnull__ (2, 3)))
inline void bmm_sort_heapsort(size_t const nmemb,
    int (*const compar)(size_t, size_t, void *),
    void (*const swap)(size_t, size_t, void *), void *cls) {
  for (size_t i = 1; i < nmemb / 2 + 1; ++i) {
    size_t p = nmemb / 2 - i;

    for ever {
      if (p >= nmemb / 2)
        break;

      size_t const cl = p * 2 + 1;
      size_t const cr = cl + 1;
      size_t const c = cr >= nmemb || compar(cl, cr, cls) >= 0 ? cl : cr;

      if (compar(p, c, cls) >= 0)
        break;

      swap(p, c, cls);

      p = c;
    }
  }

  for (size_t i = 1; i < nmemb; ++i) {
    size_t const j = nmemb - i;

    swap(0, j, cls);

    size_t p = 0;

    for ever {
      if (p >= j / 2)
        break;

      size_t cl = p * 2 + 1;
      size_t cr = cl + 1;
      size_t const c = cr >= j || compar(cl, cr, cls) >= 0 ? cl : cr;

      if (compar(p, c, cls) >= 0)
        break;

      swap(p, c, cls);

      p = c;
    }
  }
}

#endif
