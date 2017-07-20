#include <stddef.h>

#include "ext.h"

/// The call `bmm_swap(x, y)`
/// exchanges `x` with `y`.
/// Overflows are handled appropriately.
__attribute__ ((__nonnull__))
inline void inst(bmm_swap, A)(A *restrict const x, A *restrict const y) {
  A const tmp = *x;
  *x = *y;
  *y = tmp;
}

/// The call `bmm_map(nmemb, proc)`
/// maps over `nmemb` items with the procedure `proc`.
/// Overflows are handled appropriately.
__attribute__ ((__nonnull__ (2)))
inline void inst(bmm_map, A)(size_t const nmemb, void (*const proc)(size_t)) {
  for (size_t i = 0; i < nmemb; ++i)
    proc(i);
}

/// The call `bmm_map_cls(nmemb, proc, cls)`
/// maps over `nmemb` items with the procedure `proc`.
/// The closure `cls` is passed through to `proc`.
/// Overflows are handled appropriately.
__attribute__ ((__nonnull__ (2)))
inline void inst(bmm_map_cls, A)(size_t const nmemb,
    void (*const proc)(size_t, void *), void *const cls) {
  for (size_t i = 0; i < nmemb; ++i)
    proc(i, cls);
}

/// The call `bmm_foldl(nmemb, proc, init)`
/// folds over `nmemb` items with the procedure `proc`.
/// by starting from the left with the value `init`.
/// Overflows are handled appropriately.
__attribute__ ((__nonnull__ (2)))
inline A inst(bmm_foldl, A)(size_t const nmemb,
    A (*const proc)(size_t, A), A const init) {
  A x = init;

  for (size_t i = 0; i < nmemb; ++i)
    x = proc(i, x);

  return x;
}

/// The call `bmm_foldl_cls(nmemb, proc, init, cls)`
/// folds over `nmemb` items with the procedure `proc`.
/// by starting from the left with the value `init`.
/// The closure `cls` is passed through to `proc`.
/// Overflows are handled appropriately.
__attribute__ ((__nonnull__ (2)))
inline A inst(bmm_foldl_cls, A)(size_t const nmemb,
    A (*const proc)(size_t, A, void *), A const init, void *const cls) {
  A x = init;

  for (size_t i = 0; i < nmemb; ++i)
    x = proc(i, x, cls);

  return x;
}

/// The call `bmm_foldr(nmemb, proc, init)`
/// folds over `nmemb` items with the procedure `proc`.
/// by starting from the right with the value `init`.
/// Overflows are handled appropriately.
__attribute__ ((__nonnull__ (2)))
inline A inst(bmm_foldr, A)(size_t const nmemb,
    A (*const proc)(size_t, A), A const init) {
  A x = init;

  size_t const k = nmemb - 1;
  for (size_t i = 0; i < nmemb; ++i)
    x = proc(k - i, x);

  return x;
}

/// The call `bmm_foldr_cls(nmemb, proc, init, cls)`
/// folds over `nmemb` items with the procedure `proc`.
/// by starting from the right with the value `init`.
/// The closure `cls` is passed through to `proc`.
/// Overflows are handled appropriately.
__attribute__ ((__nonnull__ (2)))
inline A inst(bmm_foldr_cls, A)(size_t const nmemb,
    A (*const proc)(size_t, A, void *), A const init, void *const cls) {
  A x = init;

  size_t const k = nmemb - 1;
  for (size_t i = 0; i < nmemb; ++i)
    x = proc(k - i, x, cls);

  return x;
}
