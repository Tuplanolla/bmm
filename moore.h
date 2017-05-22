#ifndef BMM_MOORE_H
/// Moore neighborhoods with Chebyshev metrics.
#define BMM_MOORE_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "cpp.h"
#include "ext.h"
#include "size.h"

// This implementation is not as optimal as the interface allows,
// but should be good enough for most purposes.

/// The call `bmm_moore_qp(pij, i, d)`
/// checks whether the trial index `i` is in the Moore neighborhood
/// of a point in a finite
/// `d`-dimensional lattice and
/// sets `pij` to the corresponding offset.
__attribute__ ((__nonnull__))
inline bool bmm_moore_qp(size_t* const pij,
    size_t const itrial, size_t const ndim) {
  bmm_size_hc(pij, itrial, ndim, 3);

  return true;
}

/// The call `bmm_moore_q(pij, ij, i, d, n)`
/// checks whether the trial index `i` is in the Moore neighborhood
/// of the point `ij` in a finite `n`-wide
/// `d`-dimensional lattice and
/// sets `pij` to the corresponding offset.
__attribute__ ((__nonnull__))
inline bool bmm_moore_q(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const itrial,
    size_t const ndim, size_t const* restrict const nper) {
  bmm_size_hc(pij, itrial, ndim, 3);

  for (size_t idim = 0; idim < ndim; ++idim) {
    size_t const i = ij[idim] + pij[idim];

    if (i <= 0 || i > nper[idim])
      return false;
  }

  return true;
}

/// The call `bmm_moore_qcp(pij, ij, i, d, n, p)`
/// checks whether the trial index `i` is in the Moore neighborhood
/// of the point `ij` in a `p`-conditionally periodic `n`-wide
/// `d`-dimensional lattice and
/// sets `pij` to the corresponding offset.
__attribute__ ((__nonnull__))
inline bool bmm_moore_qcp(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const itrial,
    size_t const ndim, size_t const* restrict const nper,
    bool const* const per) {
  bmm_size_hc(pij, itrial, ndim, 3);

  for (size_t idim = 0; idim < ndim; ++idim)
    if (!per[idim]) {
      size_t const i = ij[idim] + pij[idim];

      if (i <= 0 || i > nper[idim])
        return false;
    }

  return true;
}

/// The call `bmm_moore_np(d)`
/// returns the number of cells in the Moore neighborhood
/// of a point in a periodic `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is less than or equal to one.
/// The result follows from $N_d = 3^d$.
__attribute__ ((__const__, __pure__))
inline size_t bmm_moore_np(size_t const ndim) {
  if (ndim == 0)
    return SIZE_MAX;

  return bmm_size_pow(3, ndim);
}

/// The call `bmm_moore_npr(d)`
/// returns the number of cells in the reduced Moore neighborhood
/// of a point in a periodic `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is equal to one.
/// The result follows from $N_d = 3^d - 1$.
__attribute__ ((__const__, __pure__))
inline size_t bmm_moore_npr(size_t const ndim) {
  if (ndim == 0)
    return SIZE_MAX;

  return bmm_size_pow(3, ndim) - 1;
}

/// The call `bmm_moore_nplh(d)`
/// returns the number of cells in the lower Moore half-neighborhood
/// of a point in a periodic `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is less than or
/// whose lexicographic position is less than or greater than
/// that of the origin's.
/// The result follows from $N_d = (3^d - 1) / 2 + 1$.
__attribute__ ((__const__, __pure__))
inline size_t bmm_moore_nplh(size_t const ndim) {
  if (ndim == 0)
    return SIZE_MAX;

  return bmm_size_pow(3, ndim) / 2 + 1;
}

/// The call `bmm_moore_nplhr(d)`
/// returns the number of cells in the reduced lower Moore half-neighborhood
/// of a point in a periodic `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is equal to one and
/// whose lexicographic position is less than
/// that of the origin's.
/// The result follows from $N_d = (3^d - 1) / 2$.
__attribute__ ((__const__, __pure__))
inline size_t bmm_moore_nplhr(size_t const ndim) {
  if (ndim == 0)
    return SIZE_MAX;

  return bmm_size_pow(3, ndim) / 2;
}

/// The call `bmm_moore_npuh(d)`
/// returns the number of cells in the upper Moore half-neighborhood
/// of a point in a periodic `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is less than or equal to one and
/// whose lexicographic position is greater than
/// that of the origin's.
/// The result follows from $N_d = (3^d - 1) / 2 + 1$.
__attribute__ ((__const__, __pure__))
inline size_t bmm_moore_npuh(size_t const ndim) {
  return bmm_moore_nplh(ndim);
}

/// The call `bmm_moore_npuhr(d)`
/// returns the number of cells in the reduced upper Moore half-neighborhood
/// of a point in a periodic `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is equal to one and
/// whose lexicographic position is greater than
/// that of the origin's.
/// The result follows from $N_d = (3^d - 1) / 2$.
__attribute__ ((__const__, __pure__))
inline size_t bmm_moore_npuhr(size_t const ndim) {
  return bmm_moore_nplhr(ndim);
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

  size_t nmoore = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);

  for (size_t itrial = 0; itrial < ntrial; ++itrial)
    if (bmm_moore_q(buf, ij, itrial, ndim, nper))
      ++nmoore;

  return nmoore;
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

  size_t nmoore = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = 0; itrial < ktrial; ++itrial)
    if (bmm_moore_q(buf, ij, itrial, ndim, nper))
      ++nmoore;

  for (size_t itrial = ktrial + 1; itrial < ntrial; ++itrial)
    if (bmm_moore_q(buf, ij, itrial, ndim, nper))
      ++nmoore;

  return nmoore;
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

  size_t nmoore = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = 0; itrial < ktrial; ++itrial)
    if (bmm_moore_q(buf, ij, itrial, ndim, nper))
      ++nmoore;

  if (bmm_moore_q(buf, ij, ktrial, ndim, nper))
    ++nmoore;

  return nmoore;
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

  size_t nmoore = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = 0; itrial < ktrial; ++itrial)
    if (bmm_moore_q(buf, ij, itrial, ndim, nper))
      ++nmoore;

  return nmoore;
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

  size_t nmoore = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  if (bmm_moore_q(buf, ij, ktrial, ndim, nper))
    ++nmoore;

  for (size_t itrial = ktrial + 1; itrial < ntrial; ++itrial)
    if (bmm_moore_q(buf, ij, itrial, ndim, nper))
      ++nmoore;

  return nmoore;
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

  size_t nmoore = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = ktrial + 1; itrial < ntrial; ++itrial)
    if (bmm_moore_q(buf, ij, itrial, ndim, nper))
      ++nmoore;

  return nmoore;
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

  size_t nmoore = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);

  for (size_t itrial = 0; itrial < ntrial; ++itrial)
    if (bmm_moore_qcp(buf, ij, itrial, ndim, nper, per))
      ++nmoore;

  return nmoore;
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

  size_t nmoore = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = 0; itrial < ktrial; ++itrial)
    if (bmm_moore_qcp(buf, ij, itrial, ndim, nper, per))
      ++nmoore;

  for (size_t itrial = ktrial + 1; itrial < ntrial; ++itrial)
    if (bmm_moore_qcp(buf, ij, itrial, ndim, nper, per))
      ++nmoore;

  return nmoore;
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

  size_t nmoore = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = 0; itrial < ktrial; ++itrial)
    if (bmm_moore_qcp(buf, ij, itrial, ndim, nper, per))
      ++nmoore;

  if (bmm_moore_qcp(buf, ij, ktrial, ndim, nper, per))
    ++nmoore;

  return nmoore;
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

  size_t nmoore = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = 0; itrial < ktrial; ++itrial)
    if (bmm_moore_qcp(buf, ij, itrial, ndim, nper, per))
      ++nmoore;

  return nmoore;
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

  size_t nmoore = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  if (bmm_moore_qcp(buf, ij, ktrial, ndim, nper, per))
    ++nmoore;

  for (size_t itrial = ktrial + 1; itrial < ntrial; ++itrial)
    if (bmm_moore_qcp(buf, ij, itrial, ndim, nper, per))
      ++nmoore;

  return nmoore;
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

  size_t nmoore = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = ktrial + 1; itrial < ntrial; ++itrial)
    if (bmm_moore_qcp(buf, ij, itrial, ndim, nper, per))
      ++nmoore;

  return nmoore;
}

/// The call `bmm_moore_ijp(pij, ij, i, d, n)`
/// sets the point `pij` to the Moore neighborhood cell `i`
/// of the point `ij` in a periodic `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_moore_np` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_moore_ijp(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const imoore,
    size_t const ndim, size_t const* restrict const nper) {
  size_t jmoore = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);

  for (size_t itrial = 0; itrial < ntrial; ++itrial)
    if (bmm_moore_qp(pij, itrial, ndim)) {
      if (jmoore == imoore)
        goto br;
      else
        ++jmoore;
    }

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    pij[idim] = bmm_size_dec(ij[idim] + pij[idim], 1, nper[idim] - 1);
}

/// The call `bmm_moore_ijpr(pij, ij, i, d, n)`
/// sets the point `pij` to the reduced Moore neighborhood cell `i`
/// of the point `ij` in a periodic `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_moore_npr` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_moore_ijpr(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const imoore,
    size_t const ndim, size_t const* restrict const nper) {
  size_t jmoore = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = 0; itrial < ktrial; ++itrial)
    if (bmm_moore_qp(pij, itrial, ndim)) {
      if (jmoore == imoore)
        goto br;
      else
        ++jmoore;
    }

  for (size_t itrial = ktrial + 1; itrial < ntrial; ++itrial)
    if (bmm_moore_qp(pij, itrial, ndim)) {
      if (jmoore == imoore)
        goto br;
      else
        ++jmoore;
    }

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    pij[idim] = bmm_size_dec(ij[idim] + pij[idim], 1, nper[idim] - 1);
}

/// The call `bmm_moore_ijplh(pij, ij, i, d, n)`
/// sets the point `pij` to the lower Moore half-neighborhood cell `i`
/// of the point `ij` in a periodic `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_moore_nplh` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_moore_ijplh(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const imoore,
    size_t const ndim, size_t const* restrict const nper) {
  size_t jmoore = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = 0; itrial < ktrial; ++itrial)
    if (bmm_moore_qp(pij, itrial, ndim)) {
      if (jmoore == imoore)
        goto br;
      else
        ++jmoore;
    }

  if (bmm_moore_qp(pij, ktrial, ndim)) {
    if (jmoore == imoore)
      goto br;
    else
      ++jmoore;
  }

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    pij[idim] = bmm_size_dec(ij[idim] + pij[idim], 1, nper[idim] - 1);
}

/// The call `bmm_moore_ijplhr(pij, ij, i, d, n)`
/// sets the point `pij` to the reduced lower Moore half-neighborhood cell `i`
/// of the point `ij` in a periodic `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_moore_nplhr` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_moore_ijplhr(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const imoore,
    size_t const ndim, size_t const* restrict const nper) {
  size_t jmoore = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = 0; itrial < ktrial; ++itrial)
    if (bmm_moore_qp(pij, itrial, ndim)) {
      if (jmoore == imoore)
        goto br;
      else
        ++jmoore;
    }

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    pij[idim] = bmm_size_dec(ij[idim] + pij[idim], 1, nper[idim] - 1);
}

/// The call `bmm_moore_ijpuh(pij, ij, i, d, n)`
/// sets the point `pij` to the upper Moore half-neighborhood cell `i`
/// of the point `ij` in a periodic `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_moore_npuh` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_moore_ijpuh(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const imoore,
    size_t const ndim, size_t const* restrict const nper) {
  size_t jmoore = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  if (bmm_moore_qp(pij, ktrial, ndim)) {
    if (jmoore == imoore)
      goto br;
    else
      ++jmoore;
  }

  for (size_t itrial = ktrial + 1; itrial < ntrial; ++itrial)
    if (bmm_moore_qp(pij, itrial, ndim)) {
      if (jmoore == imoore)
        goto br;
      else
        ++jmoore;
    }

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    pij[idim] = bmm_size_dec(ij[idim] + pij[idim], 1, nper[idim] - 1);
}

/// The call `bmm_moore_ijpuhr(pij, ij, i, d, n)`
/// sets the point `pij` to the reduced upper Moore half-neighborhood cell `i`
/// of the point `ij` in a periodic `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_moore_npuhr` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_moore_ijpuhr(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const imoore,
    size_t const ndim, size_t const* restrict const nper) {
  size_t jmoore = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = ktrial + 1; itrial < ntrial; ++itrial)
    if (bmm_moore_qp(pij, itrial, ndim)) {
      if (jmoore == imoore)
        goto br;
      else
        ++jmoore;
    }

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    pij[idim] = bmm_size_dec(ij[idim] + pij[idim], 1, nper[idim] - 1);
}

/// The call `bmm_moore_ij(pij, ij, i, d, n)`
/// sets the point `pij` to the Moore neighborhood cell `i`
/// of the point `ij` in a finite `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_moore_n` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_moore_ij(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const imoore,
    size_t const ndim, size_t const* restrict const nper) {
  size_t jmoore = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);

  for (size_t itrial = 0; itrial < ntrial; ++itrial)
    if (bmm_moore_q(pij, ij, itrial, ndim, nper)) {
      if (jmoore == imoore)
        goto br;
      else
        ++jmoore;
    }

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    pij[idim] = ij[idim] + pij[idim] - 1;
}

/// The call `bmm_moore_ijr(pij, ij, i, d, n)`
/// sets the point `pij` to the reduced Moore neighborhood cell `i`
/// of the point `ij` in a finite `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_moore_nr` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_moore_ijr(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const imoore,
    size_t const ndim, size_t const* restrict const nper) {
  size_t jmoore = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = 0; itrial < ktrial; ++itrial)
    if (bmm_moore_q(pij, ij, itrial, ndim, nper)) {
      if (jmoore == imoore)
        goto br;
      else
        ++jmoore;
    }

  for (size_t itrial = ktrial + 1; itrial < ntrial; ++itrial)
    if (bmm_moore_q(pij, ij, itrial, ndim, nper)) {
      if (jmoore == imoore)
        goto br;
      else
        ++jmoore;
    }

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    pij[idim] = ij[idim] + pij[idim] - 1;
}

/// The call `bmm_moore_ijlh(pij, ij, i, d, n)`
/// sets the point `pij` to the lower Moore half-neighborhood cell `i`
/// of the point `ij` in a finite `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_moore_nlh` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_moore_ijlh(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const imoore,
    size_t const ndim, size_t const* restrict const nper) {
  size_t jmoore = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = 0; itrial < ktrial; ++itrial)
    if (bmm_moore_q(pij, ij, itrial, ndim, nper)) {
      if (jmoore == imoore)
        goto br;
      else
        ++jmoore;
    }

  if (bmm_moore_q(pij, ij, ktrial, ndim, nper)) {
    if (jmoore == imoore)
      goto br;
    else
      ++jmoore;
  }

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    pij[idim] = ij[idim] + pij[idim] - 1;
}

/// The call `bmm_moore_ijlhr(pij, ij, i, d, n)`
/// sets the point `pij` to the reduced lower Moore half-neighborhood cell `i`
/// of the point `ij` in a finite `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_moore_nlhr` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_moore_ijlhr(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const imoore,
    size_t const ndim, size_t const* restrict const nper) {
  size_t jmoore = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = 0; itrial < ktrial; ++itrial)
    if (bmm_moore_q(pij, ij, itrial, ndim, nper)) {
      if (jmoore == imoore)
        goto br;
      else
        ++jmoore;
    }

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    pij[idim] = ij[idim] + pij[idim] - 1;
}

/// The call `bmm_moore_ijuh(pij, ij, i, d, n)`
/// sets the point `pij` to the upper Moore half-neighborhood cell `i`
/// of the point `ij` in a finite `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_moore_nuh` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_moore_ijuh(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const imoore,
    size_t const ndim, size_t const* restrict const nper) {
  size_t jmoore = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  if (bmm_moore_q(pij, ij, ktrial, ndim, nper)) {
    if (jmoore == imoore)
      goto br;
    else
      ++jmoore;
  }

  for (size_t itrial = ktrial + 1; itrial < ntrial; ++itrial)
    if (bmm_moore_q(pij, ij, itrial, ndim, nper)) {
      if (jmoore == imoore)
        goto br;
      else
        ++jmoore;
    }

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    pij[idim] = ij[idim] + pij[idim] - 1;
}

/// The call `bmm_moore_ijuhr(pij, ij, i, d, n)`
/// sets the point `pij` to the reduced upper Moore half-neighborhood cell `i`
/// of the point `ij` in a finite `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_moore_nuhr` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_moore_ijuhr(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const imoore,
    size_t const ndim, size_t const* restrict const nper) {
  size_t jmoore = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = ktrial + 1; itrial < ntrial; ++itrial)
    if (bmm_moore_q(pij, ij, itrial, ndim, nper)) {
      if (jmoore == imoore)
        goto br;
      else
        ++jmoore;
    }

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    pij[idim] = ij[idim] + pij[idim] - 1;
}

/// The call `bmm_moore_ijcp(pij, ij, i, d, n, p)`
/// sets the point `pij` to the Moore neighborhood cell `i`
/// of the point `ij` in a `p`-conditionally periodic `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_moore_ncp` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_moore_ijcp(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const imoore,
    size_t const ndim, size_t const* restrict const nper,
    bool const* const per) {
  size_t jmoore = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);

  for (size_t itrial = 0; itrial < ntrial; ++itrial)
    if (bmm_moore_qcp(pij, ij, itrial, ndim, nper, per)) {
      if (jmoore == imoore)
        goto br;
      else
        ++jmoore;
    }

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    if (per[idim])
      pij[idim] = bmm_size_dec(ij[idim] + pij[idim], 1, nper[idim] - 1);
    else
      pij[idim] = ij[idim] + pij[idim] - 1;
}

/// The call `bmm_moore_ijcpr(pij, ij, i, d, n, p)`
/// sets the point `pij` to the reduced Moore neighborhood cell `i`
/// of the point `ij` in a `p`-conditionally periodic `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_moore_ncpr` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_moore_ijcpr(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const imoore,
    size_t const ndim, size_t const* restrict const nper,
    bool const* const per) {
  size_t jmoore = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = 0; itrial < ktrial; ++itrial)
    if (bmm_moore_qcp(pij, ij, itrial, ndim, nper, per)) {
      if (jmoore == imoore)
        goto br;
      else
        ++jmoore;
    }

  for (size_t itrial = ktrial + 1; itrial < ntrial; ++itrial)
    if (bmm_moore_qcp(pij, ij, itrial, ndim, nper, per)) {
      if (jmoore == imoore)
        goto br;
      else
        ++jmoore;
    }

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    if (per[idim])
      pij[idim] = bmm_size_dec(ij[idim] + pij[idim], 1, nper[idim] - 1);
    else
      pij[idim] = ij[idim] + pij[idim] - 1;
}

/// The call `bmm_moore_ijcplh(pij, ij, i, d, n, p)`
/// sets the point `pij` to the lower Moore half-neighborhood cell `i`
/// of the point `ij` in a `p`-conditionally periodic `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_moore_ncplh` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_moore_ijcplh(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const imoore,
    size_t const ndim, size_t const* restrict const nper,
    bool const* const per) {
  size_t jmoore = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = 0; itrial < ktrial; ++itrial)
    if (bmm_moore_qcp(pij, ij, itrial, ndim, nper, per)) {
      if (jmoore == imoore)
        goto br;
      else
        ++jmoore;
    }

  if (bmm_moore_qcp(pij, ij, ktrial, ndim, nper, per)) {
    if (jmoore == imoore)
      goto br;
    else
      ++jmoore;
  }

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    if (per[idim])
      pij[idim] = bmm_size_dec(ij[idim] + pij[idim], 1, nper[idim] - 1);
    else
      pij[idim] = ij[idim] + pij[idim] - 1;
}

/// The call `bmm_moore_ijcplhr(pij, ij, i, d, n, p)`
/// sets the point `pij` to the reduced lower Moore half-neighborhood cell `i`
/// of the point `ij` in a `p`-conditionally periodic `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_moore_ncplhr` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_moore_ijcplhr(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const imoore,
    size_t const ndim, size_t const* restrict const nper,
    bool const* const per) {
  size_t jmoore = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = 0; itrial < ktrial; ++itrial)
    if (bmm_moore_qcp(pij, ij, itrial, ndim, nper, per)) {
      if (jmoore == imoore)
        goto br;
      else
        ++jmoore;
    }

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    if (per[idim])
      pij[idim] = bmm_size_dec(ij[idim] + pij[idim], 1, nper[idim] - 1);
    else
      pij[idim] = ij[idim] + pij[idim] - 1;
}

/// The call `bmm_moore_ijcpuh(pij, ij, i, d, n, p)`
/// sets the point `pij` to the upper Moore half-neighborhood cell `i`
/// of the point `ij` in a `p`-conditionally periodic `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_moore_ncpuh` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_moore_ijcpuh(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const imoore,
    size_t const ndim, size_t const* restrict const nper,
    bool const* const per) {
  size_t jmoore = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  if (bmm_moore_qcp(pij, ij, ktrial, ndim, nper, per)) {
    if (jmoore == imoore)
      goto br;
    else
      ++jmoore;
  }

  for (size_t itrial = ktrial + 1; itrial < ntrial; ++itrial)
    if (bmm_moore_qcp(pij, ij, itrial, ndim, nper, per)) {
      if (jmoore == imoore)
        goto br;
      else
        ++jmoore;
    }

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    if (per[idim])
      pij[idim] = bmm_size_dec(ij[idim] + pij[idim], 1, nper[idim] - 1);
    else
      pij[idim] = ij[idim] + pij[idim] - 1;
}

/// The call `bmm_moore_ijcpuhr(pij, ij, i, d, n, p)`
/// sets the point `pij` to the reduced upper Moore half-neighborhood cell `i`
/// of the point `ij` in a `p`-conditionally periodic `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_moore_ncpuhr` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_moore_ijcpuhr(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const imoore,
    size_t const ndim, size_t const* restrict const nper,
    bool const* const per) {
  size_t jmoore = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = ktrial + 1; itrial < ntrial; ++itrial)
    if (bmm_moore_qcp(pij, ij, itrial, ndim, nper, per)) {
      if (jmoore == imoore)
        goto br;
      else
        ++jmoore;
    }

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    if (per[idim])
      pij[idim] = bmm_size_dec(ij[idim] + pij[idim], 1, nper[idim] - 1);
    else
      pij[idim] = ij[idim] + pij[idim] - 1;
}

#endif
