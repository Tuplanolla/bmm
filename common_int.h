/// Common operations for integer types.

#include <stddef.h>

#include "ext.h"

/// This structure holds the quotient and remainder of division
/// in unspecified order.
typedef struct {
  A quot;
  A rem;
} inst(bmm_div_t, A);

/// The call `inst(bmm_div, A)(x, y)`
/// returns the quotient and remainder of `x` divided by `y`
/// in `qr` such that `y * qr.quot + qr.rem == x`,
/// Overflows are impossible both internally and externally.
__attribute__ ((__const__, __pure__))
inline inst(bmm_div_t, A) inst(bmm_div, A)(A const x, A const y) {
  inst(bmm_div_t, A) const qr = {
    .quot = x / y,
    .rem = x % y
  };

  return qr;
}

/// The call `bmm_hc(pij, i, ndim, nper)`
/// sets the index vector `pij` to the index `i`
/// in a hypercube with dimension `ndim` and side length `nper`.
/// Overflows are impossible both internally and externally.
__attribute__ ((__nonnull__))
inline void inst(bmm_hc, A)(A *const pij, A const i,
    size_t const ndim, A const nper) {
  inst(bmm_div_t, A) qr = {.quot = i};
  for (size_t idim = 0; idim < ndim; ++idim) {
    qr = inst(bmm_div, A)(qr.quot, nper);

    pij[ndim - 1 - idim] = qr.rem;
  }
}

/// The call `bmm_unhc(ij, ndim, nper)`
/// returns the index of the index vector `ij`
/// in a hypercube with dimension `ndim` and side length `nper`.
/// Overflows are impossible internally but possible externally.
__attribute__ ((__nonnull__, __pure__))
inline A inst(bmm_unhc, A)(A const *const ij,
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
inline void inst(bmm_hcd, A)(A *restrict const pij, A const i,
    size_t const ndim, A const *restrict const nper) {
  inst(bmm_div_t, A) qr = {.quot = i};
  for (size_t idim = 0; idim < ndim; ++idim) {
    qr = inst(bmm_div, A)(qr.quot, nper[ndim - 1 - idim]);

    pij[ndim - 1 - idim] = qr.rem;
  }

  // The following implementation is suitable for loop fusion,
  // but less reliable.
  // size_t *const buf = alloca(ndim * sizeof *buf);
  //
  // inst(bmm_div_t, A) qr = {.quot = i};
  // for (size_t idim = 0; idim < ndim; ++idim) {
  //   qr = inst(bmm_div, A)(qr.quot, nper[ndim - 1 - idim]);
  //
  //   buf[ndim - 1 - idim] = qr.rem;
  // }
  //
  // for (size_t idim = 0; idim < ndim; ++idim)
  //   pij[idim] = buf[idim];

  // The following implementation is suitable for loop fusion,
  // but slower (quadratic instead of linear).
  // for (size_t idim = 0; idim < ndim; ++idim) {
  //   inst(bmm_div_t, A) qr = {.quot = i};
  //   for (size_t jdim = 0; jdim < ndim - idim; ++jdim)
  //     qr = inst(bmm_div, A)(qr.quot, nper[ndim - 1 - jdim]);
  //
  //   pij[idim] = qr.rem;
  // }
}

/// The call `bmm_unhc(ij, ndim, nper)`
/// returns the index of the index vector `ij`
/// in a hypercuboid with dimension `ndim` and side lengths `nper`.
/// Overflows are impossible internally but possible externally.
__attribute__ ((__nonnull__, __pure__))
inline A inst(bmm_unhcd, A)(A const *restrict const ij,
    size_t const ndim, A const *restrict const nper) {
  A i = 0;

  for (size_t idim = 0; idim < ndim; ++idim) {
    i *= nper[idim];
    i += ij[idim];
  }

  return i;
}
