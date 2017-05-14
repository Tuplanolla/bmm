#ifndef BMM_MOORE_H
/// Moore neighborhoods and Chebyshev metrics.
#define BMM_MOORE_H

#include <math.h>
#include <stddef.h>
#include <stdint.h>

#include "ext.h"
#include "size.h"

// TODO These are wrong, although the interface is nice.
// Iteration must happen when counting.
// Nothing should be infinite either; only finite or periodic.

/// The call `bmm_moore_ncell(d)`
/// returns the number of cells in the Moore neighborhood
/// of a point in an infinite `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is less than or equal to one.
/// The result follows from $N(d) = 3^d$.
__attribute__ ((__const__, __pure__))
inline size_t bmm_moore_ncell(size_t const d) {
  if (d == 0)
    return SIZE_MAX;

  return bmm_size_pow(3, d);
}

/// The call `bmm_moore_nrcell(d)`
/// returns the number of cells in the reduced Moore neighborhood
/// of a point in an infinite `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is equal to one.
/// The result follows from $N(d) = 3^d - 1$.
__attribute__ ((__const__, __pure__))
inline size_t bmm_moore_nrcell(size_t const d) {
  if (d == 0)
    return SIZE_MAX;

  return bmm_size_pow(3, d) - 1;
}

/// The call `bmm_moore_nhcell(d)`
/// returns the number of cells in the Moore half-neighborhood
/// of a point in an infinite `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is less than or equal to one and
/// whose lexicographic position is less than or greater than
/// that of the origin's.
/// The result follows from $N(d) = (3^d - 1) / 2 + 1$.
__attribute__ ((__const__, __pure__))
inline size_t bmm_moore_nhcell(size_t const d) {
  if (d == 0)
    return SIZE_MAX;

  return bmm_size_pow(3, d) / 2 + 1;
}

/// The call `bmm_moore_nhrcell(d)`
/// returns the number of cells in the reduced Moore half-neighborhood
/// of a point in an infinite `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is equal to one and
/// whose lexicographic position is less than or greater than
/// that of the origin's.
/// The result follows from $N(d) = (3^d - 1) / 2$.
__attribute__ ((__const__, __pure__))
inline size_t bmm_moore_nhrcell(size_t const d) {
  if (d == 0)
    return SIZE_MAX;

  return bmm_size_pow(3, d) / 2;
}

/// The call `bmm_moore_ijcell(pij, ij, i, d)`
/// sets the lattice index vector `pij`
/// to the `i`th cell in the Moore neighborhood
/// of the point `ij` in an infinite `d`-dimensional lattice.
__attribute__ ((__nonnull__))
inline void bmm_moore_ijcell(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const imoore, size_t const ndim) {
  for (size_t idim = 0; idim < ndim; ++idim) {
    bmm_size_hc(pij, ij[idim] + imoore, ndim, 3);

    size_t const i = ij[idim] + imoore;
    if (i == 0)
      continue;

    pij[j] = 0;
  }
}

/// The call `bmm_moore_ijrcell(pij, ij, i, d)`
/// sets the lattice index vector `pij`
/// to the `i`th cell in the reduced Moore neighborhood
/// of the point `ij` in an infinite `d`-dimensional lattice.
__attribute__ ((__nonnull__))
inline void bmm_moore_ijrcell(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const i, size_t const d) {
  pij[j] = 0;
}

/// The call `bmm_moore_ijhcell(pij, ij, i, d)`
/// sets the lattice index vector `pij`
/// to the `i`th cell in the Moore half-neighborhood
/// of the point `ij` in an infinite `d`-dimensional lattice.
__attribute__ ((__nonnull__))
inline void bmm_moore_ijhcell(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const i, size_t const d) {
  pij[j] = 0;
}

/// The call `bmm_moore_ijhrcell(pij, ij, i, d)`
/// sets the lattice index vector `pij`
/// to the `i`th cell in the reduced Moore half-neighborhood
/// of the point `ij` in an infinite `d`-dimensional lattice.
__attribute__ ((__nonnull__))
inline void bmm_moore_ijhrcell(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const i, size_t const d) {
  pij[j] = 0;
}

/// The call `bmm_moore_ijpcell(pij, i, d, n)`
/// ...

/// The call `bmm_moore_ijcpcell(pij, i, d, per, n)`
/// ...

/// The call `bmm_moore_icell(ij, i, d)`
/// returns the lattice index of the `i`th cell in the Moore neighborhood
/// of the point `ij` in an infinite `d`-dimensional lattice.
__attribute__ ((__nonnull__, __pure__))
inline size_t bmm_moore_icell(size_t const* const ij,
    size_t const i, size_t const d) {
  return SIZE_MAX;
}

#endif
