#ifndef BMM_MOORE_H
/// Moore neighborhoods and Chebyshev metrics.
#define BMM_MOORE_H

#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "cpp.h"
#include "ext.h"
#include "size.h"

/// The call `bmm_moore_qp()`
/// checks whether an index is in the Moore neighborhood
/// of a point in a periodic lattice.
/// The result is trivially true,
/// but this function is still provided
/// to be consistent with `bmm_moore_q` and `bmm_moore_qcp`.
__attribute__ ((__const__, __pure__))
inline bool bmm_moore_qp(void) {
  return true;
}

/// The call `bmm_moore_np(d)`
/// returns the number of cells in the Moore neighborhood
/// of a point in a periodic `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is less than or equal to one.
/// The result follows from $N_d = 3^d$.
__attribute__ ((__const__, __pure__))
inline size_t bmm_moore_np(size_t const d) {
  if (d == 0)
    return SIZE_MAX;

  return bmm_size_pow(3, d);
}

/// The call `bmm_moore_npr(d)`
/// returns the number of cells in the reduced Moore neighborhood
/// of a point in a periodic `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is equal to one.
/// The result follows from $N_d = 3^d - 1$.
__attribute__ ((__const__, __pure__))
inline size_t bmm_moore_npr(size_t const d) {
  if (d == 0)
    return SIZE_MAX;

  return bmm_size_pow(3, d) - 1;
}

/// The call `bmm_moore_nph(d)`
/// returns the number of cells in the Moore half-neighborhood
/// of a point in a periodic `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is less than or equal to one and
/// whose lexicographic position is less than or greater than
/// that of the origin's.
/// The result follows from $N_d = (3^d - 1) / 2 + 1$.
__attribute__ ((__const__, __pure__))
inline size_t bmm_moore_nph(size_t const d) {
  if (d == 0)
    return SIZE_MAX;

  return bmm_size_pow(3, d) / 2 + 1;
}

/// The call `bmm_moore_nphr(d)`
/// returns the number of cells in the reduced Moore half-neighborhood
/// of a point in a periodic `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is equal to one and
/// whose lexicographic position is less than or greater than
/// that of the origin's.
/// The result follows from $N_d = (3^d - 1) / 2$.
__attribute__ ((__const__, __pure__))
inline size_t bmm_moore_nphr(size_t const d) {
  if (d == 0)
    return SIZE_MAX;

  return bmm_size_pow(3, d) / 2;
}

// It is not a mistake that procedures using a scratch space buffer
// are marked pure, because the buffer becomes undefined afterwards.

/// The call `bmm_moore_q(buf, ij, i, d, n)`
/// checks whether the index `i` is in the Moore neighborhood
/// of the point `ij` in a finite `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is less than or equal to one.
/// The array `buf` of length `d` is used as scratch space.
__attribute__ ((__nonnull__, __pure__))
inline bool bmm_moore_q(size_t* restrict const buf,
    size_t const* restrict const ij, size_t const imoore,
    size_t const ndim, size_t const* restrict const nper) {
  bmm_size_hc(buf, imoore, ndim, 3);

  for (size_t idim = 0; idim < ndim; ++idim) {
    size_t const i = ij[idim] + buf[idim];

    if (i <= 0 || i > nper[idim])
      return false;
  }

  return true;
}

/// The call `bmm_moore_n(buf, ij, d, n)`
/// returns the number of cells in the Moore neighborhood
/// of the point `ij` in a finite `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is less than or equal to one.
/// The array `buf` of length `d` is used as scratch space.
__attribute__ ((__nonnull__, __pure__))
inline size_t bmm_moore_n(size_t* restrict const buf,
    size_t const* restrict const ij, size_t const ndim,
    size_t const* restrict const nper) {
  if (ndim == 0)
    return SIZE_MAX;

  size_t n = 0;

  size_t const nmoore = bmm_size_pow(3, ndim);

  for (size_t imoore = 0; imoore < nmoore; ++imoore)
    if (bmm_moore_q(buf, ij, imoore, ndim, nper))
      ++n;

  return n;
}

/// The call `bmm_moore_nr(buf, ij, d, n)`
/// returns the number of cells in the reduced Moore neighborhood
/// of the point `ij` in a finite `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is equal to one.
/// The array `buf` of length `d` is used as scratch space.
__attribute__ ((__nonnull__, __pure__))
inline size_t bmm_moore_nr(size_t* restrict const buf,
    size_t const* restrict const ij, size_t const ndim,
    size_t const* restrict const nper) {
  if (ndim == 0)
    return SIZE_MAX;

  size_t n = 0;

  size_t const nmoore = bmm_size_pow(3, ndim);
  size_t const kmoore = nmoore / 2;

  for (size_t imoore = 0; imoore < kmoore; ++imoore)
    if (bmm_moore_q(buf, ij, imoore, ndim, nper))
      ++n;

  for (size_t imoore = kmoore + 1; imoore < nmoore; ++imoore)
    if (bmm_moore_q(buf, ij, imoore, ndim, nper))
      ++n;

  return n;
}

/// The call `bmm_moore_nlh(buf, ij, d, n)`
/// returns the number of cells in the lower Moore half-neighborhood
/// of the point `ij` in a finite `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is less than or equal to one and
/// whose lexicographic position is less than or greater than
/// that of the origin's.
/// The array `buf` of length `d` is used as scratch space.
__attribute__ ((__nonnull__, __pure__))
inline size_t bmm_moore_nlh(size_t* restrict const buf,
    size_t const* restrict const ij, size_t const ndim,
    size_t const* restrict const nper) {
  if (ndim == 0)
    return SIZE_MAX;

  size_t n = 0;

  size_t const nmoore = bmm_size_pow(3, ndim);
  size_t const kmoore = nmoore / 2;

  for (size_t imoore = 0; imoore < kmoore; ++imoore)
    if (bmm_moore_q(buf, ij, imoore, ndim, nper))
      ++n;

  if (bmm_moore_q(buf, ij, kmoore, ndim, nper))
    ++n;

  return n;
}

/// The call `bmm_moore_nlhr(buf, ij, d, n)`
/// returns the number of cells in the reduced lower Moore half-neighborhood
/// of the point `ij` in a finite `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is equal to one and
/// whose lexicographic position is less than or greater than
/// that of the origin's.
/// The array `buf` of length `d` is used as scratch space.
__attribute__ ((__nonnull__, __pure__))
inline size_t bmm_moore_nlhr(size_t* restrict const buf,
    size_t const* restrict const ij, size_t const ndim,
    size_t const* restrict const nper) {
  if (ndim == 0)
    return SIZE_MAX;

  size_t n = 0;

  size_t const nmoore = bmm_size_pow(3, ndim);
  size_t const kmoore = nmoore / 2;

  for (size_t imoore = 0; imoore < kmoore; ++imoore)
    if (bmm_moore_q(buf, ij, imoore, ndim, nper))
      ++n;

  return n;
}

/// The call `bmm_moore_nuh(buf, ij, d, n)`
/// returns the number of cells in the upper Moore half-neighborhood
/// of the point `ij` in a finite `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is less than or equal to one and
/// whose lexicographic position is less than or greater than
/// that of the origin's.
/// The array `buf` of length `d` is used as scratch space.
__attribute__ ((__nonnull__, __pure__))
inline size_t bmm_moore_nuh(size_t* restrict const buf,
    size_t const* restrict const ij, size_t const ndim,
    size_t const* restrict const nper) {
  if (ndim == 0)
    return SIZE_MAX;

  size_t n = 0;

  size_t const nmoore = bmm_size_pow(3, ndim);
  size_t const kmoore = nmoore / 2;

  if (bmm_moore_q(buf, ij, kmoore, ndim, nper))
    ++n;

  for (size_t imoore = kmoore + 1; imoore < nmoore; ++imoore)
    if (bmm_moore_q(buf, ij, imoore, ndim, nper))
      ++n;

  return n;
}

/// The call `bmm_moore_nuhr(buf, ij, d, n)`
/// returns the number of cells in the reduced upper Moore half-neighborhood
/// of the point `ij` in a finite `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is equal to one and
/// whose lexicographic position is less than or greater than
/// that of the origin's.
/// The array `buf` of length `d` is used as scratch space.
__attribute__ ((__nonnull__, __pure__))
inline size_t bmm_moore_nuhr(size_t* restrict const buf,
    size_t const* restrict const ij, size_t const ndim,
    size_t const* restrict const nper) {
  if (ndim == 0)
    return SIZE_MAX;

  size_t n = 0;

  size_t const nmoore = bmm_size_pow(3, ndim);
  size_t const kmoore = nmoore / 2;

  for (size_t imoore = kmoore + 1; imoore < nmoore; ++imoore)
    if (bmm_moore_q(buf, ij, imoore, ndim, nper))
      ++n;

  return n;
}

// TODO Variant `cp` and the rest.

/// The call `bmm_moore_ij(pij, ij, i, d)`
/// sets the lattice index vector `pij`
/// to the `i`th cell in the Moore neighborhood
/// of the point `ij` in a periodic `d`-dimensional lattice.
__attribute__ ((__nonnull__))
inline void bmm_moore_ij(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const i, size_t const ndim) {
  pij[i] = 0;
}

/// The call `bmm_moore_ijr(pij, ij, i, d)`
/// sets the lattice index vector `pij`
/// to the `i`th cell in the reduced Moore neighborhood
/// of the point `ij` in a periodic `d`-dimensional lattice.
__attribute__ ((__nonnull__))
inline void bmm_moore_ijr(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const i, size_t const d) {
  pij[i] = 0;
}

/// The call `bmm_moore_ijh(pij, ij, i, d)`
/// sets the lattice index vector `pij`
/// to the `i`th cell in the Moore half-neighborhood
/// of the point `ij` in a periodic `d`-dimensional lattice.
__attribute__ ((__nonnull__))
inline void bmm_moore_ijh(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const i, size_t const d) {
  pij[i] = 0;
}

/// The call `bmm_moore_ijhr(pij, ij, i, d)`
/// sets the lattice index vector `pij`
/// to the `i`th cell in the reduced Moore half-neighborhood
/// of the point `ij` in a periodic `d`-dimensional lattice.
__attribute__ ((__nonnull__))
inline void bmm_moore_ijhr(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const i, size_t const d) {
  pij[i] = 0;
}

/// The call `bmm_moore_ijp(pij, i, d, n)`
/// ...

/// The call `bmm_moore_ijcp(pij, i, d, per, n)`
/// ...

/// The call `bmm_moore_i(ij, i, d)`
/// returns the lattice index of the `i`th cell in the Moore neighborhood
/// of the point `ij` in a periodic `d`-dimensional lattice.
__attribute__ ((__nonnull__, __pure__))
inline size_t bmm_moore_i(size_t const* const ij,
    size_t const i, size_t const d) {
  return SIZE_MAX;
}

#endif
