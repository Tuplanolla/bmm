/// Common operations for unsigned integer types.

#include <stddef.h>

#include "ext.h"

/// The call `bmm_wrap(x, a, b)`
/// finds such `y` that `a <= y < b`
/// by shifting `x` by the appropriate number of multiples of `b - a`.
/// If `b <= a`, the behavior is undefined.
/// Overflows are impossible both internally and externally.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A type(bmm_wrap, A)(A const x, A const a, A const b) {
#ifndef DEBUG
  dynamic_assert(b > a, "Invalid argument");
#endif

  A const c = b - a;
  A const xmc = x % c;
  A const amc = a % c;

  return (xmc >= amc ? xmc - amc : c - (amc - xmc)) + a;

  // The following implementation is easier to understand,
  // but susceptible to overflowing.
  // A const c = b - a;
  //
  // return (x - a) % c + a;

  // The following implementation is easier to understand,
  // but slower (linear instead of constant).
  // A const c = b - a;
  //
  // A y = x;
  //
  // if (y < a)
  //   do
  //     y += c;
  //   while (y < a);
  // else if (y >= b)
  //   do
  //     y -= c;
  //   while (y >= b);
  //
  // return y;
}

/// The call `bmm_hc(pij, i, ndim, nper)`
/// sets the index vector `pij` to the index `i`
/// in a hypercube with dimension `ndim` and side length `nper`.
/// Overflows are impossible both internally and externally.
__attribute__ ((__nonnull__))
inline void type(bmm_hc, A)(A *const pij, A const i,
    size_t const ndim, A const nper) {
  type(bmm_div_t, A) dm = {.div = i};
  for (size_t idim = 0; idim < ndim; ++idim) {
    dm = type(bmm_div, A)(dm.div, nper);

    pij[ndim - 1 - idim] = dm.mod;
  }
}

/// The call `bmm_unhc(ij, ndim, nper)`
/// returns the index of the index vector `ij`
/// in a hypercube with dimension `ndim` and side length `nper`.
/// Overflows are impossible internally but possible externally.
__attribute__ ((__nonnull__, __pure__))
inline A type(bmm_unhc, A)(A const *const ij,
    size_t const ndim, A const nper) {
  A i = 0;

  for (size_t idim = 0; idim < ndim; ++idim) {
    i *= nper;
    i += ij[idim];
  }

  return i;
}

/// The call `bmm_hc(pij, i, ndim, nper)`
/// sets the index vector `pij` to the index `i`
/// in a hypercuboid with dimension `ndim` and side lengths `nper`.
/// Overflows are impossible both internally and externally.
__attribute__ ((__nonnull__))
inline void type(bmm_hcd, A)(A *restrict const pij, A const i,
    size_t const ndim, A const *restrict const nper) {
  type(bmm_div_t, A) dm = {.div = i};
  for (size_t idim = 0; idim < ndim; ++idim) {
    dm = type(bmm_div, A)(dm.div, nper[ndim - 1 - idim]);

    pij[ndim - 1 - idim] = dm.mod;
  }

  // The following implementation is suitable for loop fusion,
  // but less reliable.
  // size_t *const buf = alloca(ndim * sizeof *buf);
  //
  // type(bmm_div_t, A) dm = {.div = i};
  // for (size_t idim = 0; idim < ndim; ++idim) {
  //   dm = type(bmm_div, A)(dm.div, nper[ndim - 1 - idim]);
  //
  //   buf[ndim - 1 - idim] = dm.mod;
  // }
  //
  // for (size_t idim = 0; idim < ndim; ++idim)
  //   pij[idim] = buf[idim];

  // The following implementation is suitable for loop fusion,
  // but slower (quadratic instead of linear).
  // for (size_t idim = 0; idim < ndim; ++idim) {
  //   type(bmm_div_t, A) dm = {.div = i};
  //   for (size_t jdim = 0; jdim < ndim - idim; ++jdim)
  //     dm = type(bmm_div, A)(dm.div, nper[ndim - 1 - jdim]);
  //
  //   pij[idim] = dm.mod;
  // }
}

/// The call `bmm_unhc(ij, ndim, nper)`
/// returns the index of the index vector `ij`
/// in a hypercuboid with dimension `ndim` and side lengths `nper`.
/// Overflows are impossible internally but possible externally.
__attribute__ ((__nonnull__, __pure__))
inline A type(bmm_unhcd, A)(A const *restrict const ij,
    size_t const ndim, A const *restrict const nper) {
  A i = 0;

  for (size_t idim = 0; idim < ndim; ++idim) {
    i *= nper[idim];
    i += ij[idim];
  }

  return i;
}