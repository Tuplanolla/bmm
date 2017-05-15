#ifndef BMM_MOORE_H
/// Moore neighborhoods and Chebyshev metrics.
#define BMM_MOORE_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "cpp.h"
#include "ext.h"
#include "size.h"

// This implementation is not as optimal as the interface allows,
// but should be good enough for most purposes.

/// The call `bmm_moore_qp(pij, ij, i, d, n)`
/// checks whether the index `i` is in the Moore neighborhood
/// of the point `ij` in a finite `n`-wide
/// `d`-dimensional lattice and
/// sets `pij` to the corresponding offset.
__attribute__ ((__nonnull__))
inline bool bmm_moore_qp(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const imoore,
    size_t const ndim, size_t const* restrict const nper) {
  bmm_size_hc(pij, imoore, ndim, 3);

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

/// The call `bmm_moore_q(pij, ij, i, d, n)`
/// checks whether the index `i` is in the Moore neighborhood
/// of the point `ij` in a finite `n`-wide
/// `d`-dimensional lattice and
/// sets `pij` to the corresponding offset.
__attribute__ ((__nonnull__))
inline bool bmm_moore_q(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const imoore,
    size_t const ndim, size_t const* restrict const nper) {
  bmm_size_hc(pij, imoore, ndim, 3);

  for (size_t idim = 0; idim < ndim; ++idim) {
    size_t const i = ij[idim] + pij[idim];

    if (i <= 0 || i > nper[idim])
      return false;
  }

  return true;
}

/// The call `bmm_moore_n(buf, ij, d, n)`
/// returns the number of cells in the Moore neighborhood
/// of the point `ij` in a finite `n`-wide
/// `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is less than or equal to one.
/// The array `buf` of length `d` is used as scratch space.
__attribute__ ((__nonnull__))
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
/// of the point `ij` in a finite `n`-wide
/// `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is equal to one.
/// The array `buf` of length `d` is used as scratch space.
__attribute__ ((__nonnull__))
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
/// of the point `ij` in a finite `n`-wide
/// `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is less than or equal to one and
/// whose lexicographic position is less than
/// that of the origin's.
/// The array `buf` of length `d` is used as scratch space.
__attribute__ ((__nonnull__))
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
/// of the point `ij` in a finite `n`-wide
/// `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is equal to one and
/// whose lexicographic position is less than
/// that of the origin's.
/// The array `buf` of length `d` is used as scratch space.
__attribute__ ((__nonnull__))
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
/// of the point `ij` in a finite `n`-wide
/// `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is less than or equal to one and
/// whose lexicographic position is greater than
/// that of the origin's.
/// The array `buf` of length `d` is used as scratch space.
__attribute__ ((__nonnull__))
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
/// of the point `ij` in a finite `n`-wide
/// `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is equal to one and
/// whose lexicographic position is greater than
/// that of the origin's.
/// The array `buf` of length `d` is used as scratch space.
__attribute__ ((__nonnull__))
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

/// The call `bmm_moore_qcp(pij, ij, i, d, n, p)`
/// checks whether the index `i` is in the Moore neighborhood
/// of the point `ij` in a `p`-conditionally periodic `n`-wide
/// `d`-dimensional lattice and
/// sets `pij` to the corresponding offset.
__attribute__ ((__nonnull__))
inline bool bmm_moore_qcp(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const imoore,
    size_t const ndim, size_t const* restrict const nper,
    bool const* const per) {
  bmm_size_hc(pij, imoore, ndim, 3);

  for (size_t idim = 0; idim < ndim; ++idim)
    if (!per[idim]) {
      size_t const i = ij[idim] + pij[idim];

      if (i <= 0 || i > nper[idim])
        return false;
    }

  return true;
}

/// The call `bmm_moore_ncp(buf, ij, d, n, p)`
/// returns the number of cells in the Moore neighborhood
/// of the point `ij` in a `p`-conditionally periodic `n`-wide
/// `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is less than or equal to one.
/// The array `buf` of length `d` is used as scratch space.
__attribute__ ((__nonnull__))
inline size_t bmm_moore_ncp(size_t* restrict const buf,
    size_t const* restrict const ij, size_t const ndim,
    size_t const* restrict const nper, bool const* const per) {
  if (ndim == 0)
    return SIZE_MAX;

  size_t n = 0;

  size_t const nmoore = bmm_size_pow(3, ndim);

  for (size_t imoore = 0; imoore < nmoore; ++imoore)
    if (bmm_moore_qcp(buf, ij, imoore, ndim, nper, per))
      ++n;

  return n;
}

/// The call `bmm_moore_ncpr(buf, ij, d, n, p)`
/// returns the number of cells in the reduced Moore neighborhood
/// of the point `ij` in a `p`-conditionally periodic `n`-wide
/// `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is equal to one.
/// The array `buf` of length `d` is used as scratch space.
__attribute__ ((__nonnull__))
inline size_t bmm_moore_ncpr(size_t* restrict const buf,
    size_t const* restrict const ij, size_t const ndim,
    size_t const* restrict const nper, bool const* const per) {
  if (ndim == 0)
    return SIZE_MAX;

  size_t n = 0;

  size_t const nmoore = bmm_size_pow(3, ndim);
  size_t const kmoore = nmoore / 2;

  for (size_t imoore = 0; imoore < kmoore; ++imoore)
    if (bmm_moore_qcp(buf, ij, imoore, ndim, nper, per))
      ++n;

  for (size_t imoore = kmoore + 1; imoore < nmoore; ++imoore)
    if (bmm_moore_qcp(buf, ij, imoore, ndim, nper, per))
      ++n;

  return n;
}

/// The call `bmm_moore_ncplh(buf, ij, d, n, p)`
/// returns the number of cells in the lower Moore half-neighborhood
/// of the point `ij` in a `p`-conditionally periodic `n`-wide
/// `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is less than or equal to one and
/// whose lexicographic position is less than
/// that of the origin's.
/// The array `buf` of length `d` is used as scratch space.
__attribute__ ((__nonnull__))
inline size_t bmm_moore_ncplh(size_t* restrict const buf,
    size_t const* restrict const ij, size_t const ndim,
    size_t const* restrict const nper, bool const* const per) {
  if (ndim == 0)
    return SIZE_MAX;

  size_t n = 0;

  size_t const nmoore = bmm_size_pow(3, ndim);
  size_t const kmoore = nmoore / 2;

  for (size_t imoore = 0; imoore < kmoore; ++imoore)
    if (bmm_moore_qcp(buf, ij, imoore, ndim, nper, per))
      ++n;

  if (bmm_moore_qcp(buf, ij, kmoore, ndim, nper, per))
    ++n;

  return n;
}

/// The call `bmm_moore_ncplhr(buf, ij, d, n, p)`
/// returns the number of cells in the reduced lower Moore half-neighborhood
/// of the point `ij` in a `p`-conditionally periodic `n`-wide
/// `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is equal to one and
/// whose lexicographic position is less than
/// that of the origin's.
/// The array `buf` of length `d` is used as scratch space.
__attribute__ ((__nonnull__))
inline size_t bmm_moore_ncplhr(size_t* restrict const buf,
    size_t const* restrict const ij, size_t const ndim,
    size_t const* restrict const nper, bool const* const per) {
  if (ndim == 0)
    return SIZE_MAX;

  size_t n = 0;

  size_t const nmoore = bmm_size_pow(3, ndim);
  size_t const kmoore = nmoore / 2;

  for (size_t imoore = 0; imoore < kmoore; ++imoore)
    if (bmm_moore_qcp(buf, ij, imoore, ndim, nper, per))
      ++n;

  return n;
}

/// The call `bmm_moore_ncpuh(buf, ij, d, n, p)`
/// returns the number of cells in the upper Moore half-neighborhood
/// of the point `ij` in a `p`-conditionally periodic `n`-wide
/// `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is less than or equal to one and
/// whose lexicographic position is greater than
/// that of the origin's.
/// The array `buf` of length `d` is used as scratch space.
__attribute__ ((__nonnull__))
inline size_t bmm_moore_ncpuh(size_t* restrict const buf,
    size_t const* restrict const ij, size_t const ndim,
    size_t const* restrict const nper, bool const* const per) {
  if (ndim == 0)
    return SIZE_MAX;

  size_t n = 0;

  size_t const nmoore = bmm_size_pow(3, ndim);
  size_t const kmoore = nmoore / 2;

  if (bmm_moore_qcp(buf, ij, kmoore, ndim, nper, per))
    ++n;

  for (size_t imoore = kmoore + 1; imoore < nmoore; ++imoore)
    if (bmm_moore_qcp(buf, ij, imoore, ndim, nper, per))
      ++n;

  return n;
}

/// The call `bmm_moore_ncpuhr(buf, ij, d, n, p)`
/// returns the number of cells in the reduced upper Moore half-neighborhood
/// of the point `ij` in a `p`-conditionally periodic `n`-wide
/// `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is equal to one and
/// whose lexicographic position is greater than
/// that of the origin's.
/// The array `buf` of length `d` is used as scratch space.
__attribute__ ((__nonnull__))
inline size_t bmm_moore_ncpuhr(size_t* restrict const buf,
    size_t const* restrict const ij, size_t const ndim,
    size_t const* restrict const nper, bool const* const per) {
  if (ndim == 0)
    return SIZE_MAX;

  size_t n = 0;

  size_t const nmoore = bmm_size_pow(3, ndim);
  size_t const kmoore = nmoore / 2;

  for (size_t imoore = kmoore + 1; imoore < nmoore; ++imoore)
    if (bmm_moore_qcp(buf, ij, imoore, ndim, nper, per))
      ++n;

  return n;
}

/// The call `bmm_moore_ijp(pij, ij, i, d, n)`
/// sets the point `pij` to the index `i`
/// of the point `ij` in a periodic `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_moore_n` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_moore_ijp(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const i,
    size_t const ndim, size_t const* restrict const nper) {
  size_t j = 0;

  size_t const nmoore = bmm_size_pow(3, ndim);

  for (size_t imoore = 0; imoore < nmoore; ++imoore)
    if (bmm_moore_qp(pij, ij, imoore, ndim, nper)) {
      if (j == i)
        break;
      else
        ++j;
    }

  for (size_t idim = 0; idim < ndim; ++idim)
    pij[idim] = bmm_size_dec(ij[idim] + pij[idim], nper[idim]);
}

/// The call `bmm_moore_ij(pij, ij, i, d, n)`
/// sets the point `pij` to the index `i`
/// of the point `ij` in a finite `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_moore_n` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_moore_ij(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const i,
    size_t const ndim, size_t const* restrict const nper) {
  size_t j = 0;

  size_t const nmoore = bmm_size_pow(3, ndim);

  for (size_t imoore = 0; imoore < nmoore; ++imoore)
    if (bmm_moore_q(pij, ij, imoore, ndim, nper)) {
      if (j == i)
        break;
      else
        ++j;
    }

  for (size_t idim = 0; idim < ndim; ++idim)
    pij[idim] = ij[idim] + pij[idim] - 1;
}

#endif
