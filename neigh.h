#ifndef BMM_NEIGH_H
/// Moore neighborhoods with Chebyshev metrics.
#define BMM_NEIGH_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "cpp.h"
#include "ext.h"
#include "size.h"

// This implementation is not as optimal as the interface allows,
// but should be good enough for most purposes.

/// The call `bmm_neigh_qijp(pij, i, d)`
/// checks whether the trial index `i` is in the Moore neighborhood
/// of a point in a finite
/// `d`-dimensional lattice and
/// sets `pij` to the corresponding offset.
__attribute__ ((__nonnull__))
inline bool bmm_neigh_qijp(size_t* const pij,
    size_t const itrial, size_t const ndim) {
  bmm_size_hc(pij, itrial, ndim, 3);

  return true;
}

/// The call `bmm_neigh_qij(pij, ij, i, d, n)`
/// checks whether the trial index `i` is in the Moore neighborhood
/// of the point `ij` in a finite `n`-wide
/// `d`-dimensional lattice and
/// sets `pij` to the corresponding offset.
__attribute__ ((__nonnull__))
inline bool bmm_neigh_qij(size_t* restrict const pij,
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

/// The call `bmm_neigh_qijcp(pij, ij, i, d, n, p)`
/// checks whether the trial index `i` is in the Moore neighborhood
/// of the point `ij` in a `p`-conditionally periodic `n`-wide
/// `d`-dimensional lattice and
/// sets `pij` to the corresponding offset.
__attribute__ ((__nonnull__))
inline bool bmm_neigh_qijcp(size_t* restrict const pij,
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

/// The call `bmm_neigh_np(d)`
/// returns the number of cells in the Moore neighborhood
/// of a point in a periodic `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is less than or equal to one.
/// The result follows from $N_d = 3^d$.
__attribute__ ((__const__, __pure__))
inline size_t bmm_neigh_np(size_t const ndim) {
  if (ndim == 0)
    return SIZE_MAX;

  return bmm_size_pow(3, ndim);
}

/// The call `bmm_neigh_npr(d)`
/// returns the number of cells in the reduced Moore neighborhood
/// of a point in a periodic `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is equal to one.
/// The result follows from $N_d = 3^d - 1$.
__attribute__ ((__const__, __pure__))
inline size_t bmm_neigh_npr(size_t const ndim) {
  if (ndim == 0)
    return SIZE_MAX;

  return bmm_size_pow(3, ndim) - 1;
}

/// The call `bmm_neigh_nplh(d)`
/// returns the number of cells in the lower Moore half-neighborhood
/// of a point in a periodic `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is less than or
/// whose lexicographic position is less than or greater than
/// that of the origin's.
/// The result follows from $N_d = (3^d - 1) / 2 + 1$.
__attribute__ ((__const__, __pure__))
inline size_t bmm_neigh_nplh(size_t const ndim) {
  if (ndim == 0)
    return SIZE_MAX;

  return bmm_size_pow(3, ndim) / 2 + 1;
}

/// The call `bmm_neigh_nplhr(d)`
/// returns the number of cells in the reduced lower Moore half-neighborhood
/// of a point in a periodic `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is equal to one and
/// whose lexicographic position is less than
/// that of the origin's.
/// The result follows from $N_d = (3^d - 1) / 2$.
__attribute__ ((__const__, __pure__))
inline size_t bmm_neigh_nplhr(size_t const ndim) {
  if (ndim == 0)
    return SIZE_MAX;

  return bmm_size_pow(3, ndim) / 2;
}

/// The call `bmm_neigh_npuh(d)`
/// returns the number of cells in the upper Moore half-neighborhood
/// of a point in a periodic `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is less than or equal to one and
/// whose lexicographic position is greater than
/// that of the origin's.
/// The result follows from $N_d = (3^d - 1) / 2 + 1$.
__attribute__ ((__const__, __pure__))
inline size_t bmm_neigh_npuh(size_t const ndim) {
  return bmm_neigh_nplh(ndim);
}

/// The call `bmm_neigh_npuhr(d)`
/// returns the number of cells in the reduced upper Moore half-neighborhood
/// of a point in a periodic `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is equal to one and
/// whose lexicographic position is greater than
/// that of the origin's.
/// The result follows from $N_d = (3^d - 1) / 2$.
__attribute__ ((__const__, __pure__))
inline size_t bmm_neigh_npuhr(size_t const ndim) {
  return bmm_neigh_nplhr(ndim);
}

/// The call `bmm_neigh_n(buf, ij, d, n)`
/// returns the number of cells in the Moore neighborhood
/// of the point `ij` in a finite `n`-wide
/// `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is less than or equal to one.
/// The array `buf` of length `d` is used as scratch space.
__attribute__ ((__nonnull__))
inline size_t bmm_neigh_n(size_t* restrict const buf,
    size_t const* restrict const ij, size_t const ndim,
    size_t const* restrict const nper) {
  if (ndim == 0)
    return SIZE_MAX;

  size_t nneigh = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);

  for (size_t itrial = 0; itrial < ntrial; ++itrial)
    if (bmm_neigh_qij(buf, ij, itrial, ndim, nper))
      ++nneigh;

  return nneigh;
}

/// The call `bmm_neigh_nr(buf, ij, d, n)`
/// returns the number of cells in the reduced Moore neighborhood
/// of the point `ij` in a finite `n`-wide
/// `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is equal to one.
/// The array `buf` of length `d` is used as scratch space.
__attribute__ ((__nonnull__))
inline size_t bmm_neigh_nr(size_t* restrict const buf,
    size_t const* restrict const ij, size_t const ndim,
    size_t const* restrict const nper) {
  if (ndim == 0)
    return SIZE_MAX;

  size_t nneigh = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = 0; itrial < ktrial; ++itrial)
    if (bmm_neigh_qij(buf, ij, itrial, ndim, nper))
      ++nneigh;

  for (size_t itrial = ktrial + 1; itrial < ntrial; ++itrial)
    if (bmm_neigh_qij(buf, ij, itrial, ndim, nper))
      ++nneigh;

  return nneigh;
}

/// The call `bmm_neigh_nlh(buf, ij, d, n)`
/// returns the number of cells in the lower Moore half-neighborhood
/// of the point `ij` in a finite `n`-wide
/// `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is less than or equal to one and
/// whose lexicographic position is less than
/// that of the origin's.
/// The array `buf` of length `d` is used as scratch space.
__attribute__ ((__nonnull__))
inline size_t bmm_neigh_nlh(size_t* restrict const buf,
    size_t const* restrict const ij, size_t const ndim,
    size_t const* restrict const nper) {
  if (ndim == 0)
    return SIZE_MAX;

  size_t nneigh = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = 0; itrial < ktrial; ++itrial)
    if (bmm_neigh_qij(buf, ij, itrial, ndim, nper))
      ++nneigh;

  if (bmm_neigh_qij(buf, ij, ktrial, ndim, nper))
    ++nneigh;

  return nneigh;
}

/// The call `bmm_neigh_nlhr(buf, ij, d, n)`
/// returns the number of cells in the reduced lower Moore half-neighborhood
/// of the point `ij` in a finite `n`-wide
/// `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is equal to one and
/// whose lexicographic position is less than
/// that of the origin's.
/// The array `buf` of length `d` is used as scratch space.
__attribute__ ((__nonnull__))
inline size_t bmm_neigh_nlhr(size_t* restrict const buf,
    size_t const* restrict const ij, size_t const ndim,
    size_t const* restrict const nper) {
  if (ndim == 0)
    return SIZE_MAX;

  size_t nneigh = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = 0; itrial < ktrial; ++itrial)
    if (bmm_neigh_qij(buf, ij, itrial, ndim, nper))
      ++nneigh;

  return nneigh;
}

/// The call `bmm_neigh_nuh(buf, ij, d, n)`
/// returns the number of cells in the upper Moore half-neighborhood
/// of the point `ij` in a finite `n`-wide
/// `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is less than or equal to one and
/// whose lexicographic position is greater than
/// that of the origin's.
/// The array `buf` of length `d` is used as scratch space.
__attribute__ ((__nonnull__))
inline size_t bmm_neigh_nuh(size_t* restrict const buf,
    size_t const* restrict const ij, size_t const ndim,
    size_t const* restrict const nper) {
  if (ndim == 0)
    return SIZE_MAX;

  size_t nneigh = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  if (bmm_neigh_qij(buf, ij, ktrial, ndim, nper))
    ++nneigh;

  for (size_t itrial = ktrial + 1; itrial < ntrial; ++itrial)
    if (bmm_neigh_qij(buf, ij, itrial, ndim, nper))
      ++nneigh;

  return nneigh;
}

/// The call `bmm_neigh_nuhr(buf, ij, d, n)`
/// returns the number of cells in the reduced upper Moore half-neighborhood
/// of the point `ij` in a finite `n`-wide
/// `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is equal to one and
/// whose lexicographic position is greater than
/// that of the origin's.
/// The array `buf` of length `d` is used as scratch space.
__attribute__ ((__nonnull__))
inline size_t bmm_neigh_nuhr(size_t* restrict const buf,
    size_t const* restrict const ij, size_t const ndim,
    size_t const* restrict const nper) {
  if (ndim == 0)
    return SIZE_MAX;

  size_t nneigh = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = ktrial + 1; itrial < ntrial; ++itrial)
    if (bmm_neigh_qij(buf, ij, itrial, ndim, nper))
      ++nneigh;

  return nneigh;
}

/// The call `bmm_neigh_ncp(buf, ij, d, n, p)`
/// returns the number of cells in the Moore neighborhood
/// of the point `ij` in a `p`-conditionally periodic `n`-wide
/// `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is less than or equal to one.
/// The array `buf` of length `d` is used as scratch space.
__attribute__ ((__nonnull__))
inline size_t bmm_neigh_ncp(size_t* restrict const buf,
    size_t const* restrict const ij, size_t const ndim,
    size_t const* restrict const nper, bool const* const per) {
  if (ndim == 0)
    return SIZE_MAX;

  size_t nneigh = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);

  for (size_t itrial = 0; itrial < ntrial; ++itrial)
    if (bmm_neigh_qijcp(buf, ij, itrial, ndim, nper, per))
      ++nneigh;

  return nneigh;
}

/// The call `bmm_neigh_ncpr(buf, ij, d, n, p)`
/// returns the number of cells in the reduced Moore neighborhood
/// of the point `ij` in a `p`-conditionally periodic `n`-wide
/// `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is equal to one.
/// The array `buf` of length `d` is used as scratch space.
__attribute__ ((__nonnull__))
inline size_t bmm_neigh_ncpr(size_t* restrict const buf,
    size_t const* restrict const ij, size_t const ndim,
    size_t const* restrict const nper, bool const* const per) {
  if (ndim == 0)
    return SIZE_MAX;

  size_t nneigh = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = 0; itrial < ktrial; ++itrial)
    if (bmm_neigh_qijcp(buf, ij, itrial, ndim, nper, per))
      ++nneigh;

  for (size_t itrial = ktrial + 1; itrial < ntrial; ++itrial)
    if (bmm_neigh_qijcp(buf, ij, itrial, ndim, nper, per))
      ++nneigh;

  return nneigh;
}

/// The call `bmm_neigh_ncplh(buf, ij, d, n, p)`
/// returns the number of cells in the lower Moore half-neighborhood
/// of the point `ij` in a `p`-conditionally periodic `n`-wide
/// `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is less than or equal to one and
/// whose lexicographic position is less than
/// that of the origin's.
/// The array `buf` of length `d` is used as scratch space.
__attribute__ ((__nonnull__))
inline size_t bmm_neigh_ncplh(size_t* restrict const buf,
    size_t const* restrict const ij, size_t const ndim,
    size_t const* restrict const nper, bool const* const per) {
  if (ndim == 0)
    return SIZE_MAX;

  size_t nneigh = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = 0; itrial < ktrial; ++itrial)
    if (bmm_neigh_qijcp(buf, ij, itrial, ndim, nper, per))
      ++nneigh;

  if (bmm_neigh_qijcp(buf, ij, ktrial, ndim, nper, per))
    ++nneigh;

  return nneigh;
}

/// The call `bmm_neigh_ncplhr(buf, ij, d, n, p)`
/// returns the number of cells in the reduced lower Moore half-neighborhood
/// of the point `ij` in a `p`-conditionally periodic `n`-wide
/// `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is equal to one and
/// whose lexicographic position is less than
/// that of the origin's.
/// The array `buf` of length `d` is used as scratch space.
__attribute__ ((__nonnull__))
inline size_t bmm_neigh_ncplhr(size_t* restrict const buf,
    size_t const* restrict const ij, size_t const ndim,
    size_t const* restrict const nper, bool const* const per) {
  if (ndim == 0)
    return SIZE_MAX;

  size_t nneigh = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = 0; itrial < ktrial; ++itrial)
    if (bmm_neigh_qijcp(buf, ij, itrial, ndim, nper, per))
      ++nneigh;

  return nneigh;
}

/// The call `bmm_neigh_ncpuh(buf, ij, d, n, p)`
/// returns the number of cells in the upper Moore half-neighborhood
/// of the point `ij` in a `p`-conditionally periodic `n`-wide
/// `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is less than or equal to one and
/// whose lexicographic position is greater than
/// that of the origin's.
/// The array `buf` of length `d` is used as scratch space.
__attribute__ ((__nonnull__))
inline size_t bmm_neigh_ncpuh(size_t* restrict const buf,
    size_t const* restrict const ij, size_t const ndim,
    size_t const* restrict const nper, bool const* const per) {
  if (ndim == 0)
    return SIZE_MAX;

  size_t nneigh = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  if (bmm_neigh_qijcp(buf, ij, ktrial, ndim, nper, per))
    ++nneigh;

  for (size_t itrial = ktrial + 1; itrial < ntrial; ++itrial)
    if (bmm_neigh_qijcp(buf, ij, itrial, ndim, nper, per))
      ++nneigh;

  return nneigh;
}

/// The call `bmm_neigh_ncpuhr(buf, ij, d, n, p)`
/// returns the number of cells in the reduced upper Moore half-neighborhood
/// of the point `ij` in a `p`-conditionally periodic `n`-wide
/// `d`-dimensional lattice.
/// This includes all the cells
/// whose Chebyshev distance is equal to one and
/// whose lexicographic position is greater than
/// that of the origin's.
/// The array `buf` of length `d` is used as scratch space.
__attribute__ ((__nonnull__))
inline size_t bmm_neigh_ncpuhr(size_t* restrict const buf,
    size_t const* restrict const ij, size_t const ndim,
    size_t const* restrict const nper, bool const* const per) {
  if (ndim == 0)
    return SIZE_MAX;

  size_t nneigh = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = ktrial + 1; itrial < ntrial; ++itrial)
    if (bmm_neigh_qijcp(buf, ij, itrial, ndim, nper, per))
      ++nneigh;

  return nneigh;
}

/// The call `bmm_neigh_ijp(pij, ij, i, d, n)`
/// sets the point `pij` to the Moore neighborhood cell `i`
/// of the point `ij` in a periodic `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_neigh_np` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_neigh_ijp(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const ineigh,
    size_t const ndim, size_t const* restrict const nper) {
  size_t jneigh = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);

  for (size_t itrial = 0; itrial < ntrial; ++itrial)
    if (bmm_neigh_qijp(pij, itrial, ndim)) {
      if (jneigh == ineigh)
        goto br;
      else
        ++jneigh;
    }

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    pij[idim] = bmm_size_dec(ij[idim] + pij[idim], 1, nper[idim] - 1);
}

/// The call `bmm_neigh_ijpr(pij, ij, i, d, n)`
/// sets the point `pij` to the reduced Moore neighborhood cell `i`
/// of the point `ij` in a periodic `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_neigh_npr` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_neigh_ijpr(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const ineigh,
    size_t const ndim, size_t const* restrict const nper) {
  size_t jneigh = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = 0; itrial < ktrial; ++itrial)
    if (bmm_neigh_qijp(pij, itrial, ndim)) {
      if (jneigh == ineigh)
        goto br;
      else
        ++jneigh;
    }

  for (size_t itrial = ktrial + 1; itrial < ntrial; ++itrial)
    if (bmm_neigh_qijp(pij, itrial, ndim)) {
      if (jneigh == ineigh)
        goto br;
      else
        ++jneigh;
    }

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    pij[idim] = bmm_size_dec(ij[idim] + pij[idim], 1, nper[idim] - 1);
}

/// The call `bmm_neigh_ijplh(pij, ij, i, d, n)`
/// sets the point `pij` to the lower Moore half-neighborhood cell `i`
/// of the point `ij` in a periodic `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_neigh_nplh` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_neigh_ijplh(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const ineigh,
    size_t const ndim, size_t const* restrict const nper) {
  size_t jneigh = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = 0; itrial < ktrial; ++itrial)
    if (bmm_neigh_qijp(pij, itrial, ndim)) {
      if (jneigh == ineigh)
        goto br;
      else
        ++jneigh;
    }

  if (bmm_neigh_qijp(pij, ktrial, ndim)) {
    if (jneigh == ineigh)
      goto br;
    else
      ++jneigh;
  }

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    pij[idim] = bmm_size_dec(ij[idim] + pij[idim], 1, nper[idim] - 1);
}

/// The call `bmm_neigh_ijplhr(pij, ij, i, d, n)`
/// sets the point `pij` to the reduced lower Moore half-neighborhood cell `i`
/// of the point `ij` in a periodic `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_neigh_nplhr` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_neigh_ijplhr(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const ineigh,
    size_t const ndim, size_t const* restrict const nper) {
  size_t jneigh = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = 0; itrial < ktrial; ++itrial)
    if (bmm_neigh_qijp(pij, itrial, ndim)) {
      if (jneigh == ineigh)
        goto br;
      else
        ++jneigh;
    }

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    pij[idim] = bmm_size_dec(ij[idim] + pij[idim], 1, nper[idim] - 1);
}

/// The call `bmm_neigh_ijpuh(pij, ij, i, d, n)`
/// sets the point `pij` to the upper Moore half-neighborhood cell `i`
/// of the point `ij` in a periodic `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_neigh_npuh` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_neigh_ijpuh(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const ineigh,
    size_t const ndim, size_t const* restrict const nper) {
  size_t jneigh = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  if (bmm_neigh_qijp(pij, ktrial, ndim)) {
    if (jneigh == ineigh)
      goto br;
    else
      ++jneigh;
  }

  for (size_t itrial = ktrial + 1; itrial < ntrial; ++itrial)
    if (bmm_neigh_qijp(pij, itrial, ndim)) {
      if (jneigh == ineigh)
        goto br;
      else
        ++jneigh;
    }

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    pij[idim] = bmm_size_dec(ij[idim] + pij[idim], 1, nper[idim] - 1);
}

/// The call `bmm_neigh_ijpuhr(pij, ij, i, d, n)`
/// sets the point `pij` to the reduced upper Moore half-neighborhood cell `i`
/// of the point `ij` in a periodic `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_neigh_npuhr` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_neigh_ijpuhr(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const ineigh,
    size_t const ndim, size_t const* restrict const nper) {
  size_t jneigh = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = ktrial + 1; itrial < ntrial; ++itrial)
    if (bmm_neigh_qijp(pij, itrial, ndim)) {
      if (jneigh == ineigh)
        goto br;
      else
        ++jneigh;
    }

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    pij[idim] = bmm_size_dec(ij[idim] + pij[idim], 1, nper[idim] - 1);
}

/// The call `bmm_neigh_ij(pij, ij, i, d, n)`
/// sets the point `pij` to the Moore neighborhood cell `i`
/// of the point `ij` in a finite `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_neigh_n` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_neigh_ij(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const ineigh,
    size_t const ndim, size_t const* restrict const nper) {
  size_t jneigh = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);

  for (size_t itrial = 0; itrial < ntrial; ++itrial)
    if (bmm_neigh_qij(pij, ij, itrial, ndim, nper)) {
      if (jneigh == ineigh)
        goto br;
      else
        ++jneigh;
    }

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    pij[idim] = ij[idim] + pij[idim] - 1;
}

/// The call `bmm_neigh_ijr(pij, ij, i, d, n)`
/// sets the point `pij` to the reduced Moore neighborhood cell `i`
/// of the point `ij` in a finite `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_neigh_nr` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_neigh_ijr(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const ineigh,
    size_t const ndim, size_t const* restrict const nper) {
  size_t jneigh = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = 0; itrial < ktrial; ++itrial)
    if (bmm_neigh_qij(pij, ij, itrial, ndim, nper)) {
      if (jneigh == ineigh)
        goto br;
      else
        ++jneigh;
    }

  for (size_t itrial = ktrial + 1; itrial < ntrial; ++itrial)
    if (bmm_neigh_qij(pij, ij, itrial, ndim, nper)) {
      if (jneigh == ineigh)
        goto br;
      else
        ++jneigh;
    }

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    pij[idim] = ij[idim] + pij[idim] - 1;
}

/// The call `bmm_neigh_ijlh(pij, ij, i, d, n)`
/// sets the point `pij` to the lower Moore half-neighborhood cell `i`
/// of the point `ij` in a finite `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_neigh_nlh` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_neigh_ijlh(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const ineigh,
    size_t const ndim, size_t const* restrict const nper) {
  size_t jneigh = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = 0; itrial < ktrial; ++itrial)
    if (bmm_neigh_qij(pij, ij, itrial, ndim, nper)) {
      if (jneigh == ineigh)
        goto br;
      else
        ++jneigh;
    }

  if (bmm_neigh_qij(pij, ij, ktrial, ndim, nper)) {
    if (jneigh == ineigh)
      goto br;
    else
      ++jneigh;
  }

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    pij[idim] = ij[idim] + pij[idim] - 1;
}

/// The call `bmm_neigh_ijlhr(pij, ij, i, d, n)`
/// sets the point `pij` to the reduced lower Moore half-neighborhood cell `i`
/// of the point `ij` in a finite `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_neigh_nlhr` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_neigh_ijlhr(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const ineigh,
    size_t const ndim, size_t const* restrict const nper) {
  size_t jneigh = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = 0; itrial < ktrial; ++itrial)
    if (bmm_neigh_qij(pij, ij, itrial, ndim, nper)) {
      if (jneigh == ineigh)
        goto br;
      else
        ++jneigh;
    }

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    pij[idim] = ij[idim] + pij[idim] - 1;
}

/// The call `bmm_neigh_ijuh(pij, ij, i, d, n)`
/// sets the point `pij` to the upper Moore half-neighborhood cell `i`
/// of the point `ij` in a finite `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_neigh_nuh` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_neigh_ijuh(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const ineigh,
    size_t const ndim, size_t const* restrict const nper) {
  size_t jneigh = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  if (bmm_neigh_qij(pij, ij, ktrial, ndim, nper)) {
    if (jneigh == ineigh)
      goto br;
    else
      ++jneigh;
  }

  for (size_t itrial = ktrial + 1; itrial < ntrial; ++itrial)
    if (bmm_neigh_qij(pij, ij, itrial, ndim, nper)) {
      if (jneigh == ineigh)
        goto br;
      else
        ++jneigh;
    }

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    pij[idim] = ij[idim] + pij[idim] - 1;
}

/// The call `bmm_neigh_ijuhr(pij, ij, i, d, n)`
/// sets the point `pij` to the reduced upper Moore half-neighborhood cell `i`
/// of the point `ij` in a finite `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_neigh_nuhr` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_neigh_ijuhr(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const ineigh,
    size_t const ndim, size_t const* restrict const nper) {
  size_t jneigh = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = ktrial + 1; itrial < ntrial; ++itrial)
    if (bmm_neigh_qij(pij, ij, itrial, ndim, nper)) {
      if (jneigh == ineigh)
        goto br;
      else
        ++jneigh;
    }

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    pij[idim] = ij[idim] + pij[idim] - 1;
}

/// The call `bmm_neigh_ijcp(pij, ij, i, d, n, p)`
/// sets the point `pij` to the Moore neighborhood cell `i`
/// of the point `ij` in a `p`-conditionally periodic `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_neigh_ncp` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_neigh_ijcp(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const ineigh,
    size_t const ndim, size_t const* restrict const nper,
    bool const* const per) {
  size_t jneigh = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);

  for (size_t itrial = 0; itrial < ntrial; ++itrial)
    if (bmm_neigh_qijcp(pij, ij, itrial, ndim, nper, per)) {
      if (jneigh == ineigh)
        goto br;
      else
        ++jneigh;
    }

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    if (per[idim])
      pij[idim] = bmm_size_dec(ij[idim] + pij[idim], 1, nper[idim] - 1);
    else
      pij[idim] = ij[idim] + pij[idim] - 1;
}

/// The call `bmm_neigh_ijcpr(pij, ij, i, d, n, p)`
/// sets the point `pij` to the reduced Moore neighborhood cell `i`
/// of the point `ij` in a `p`-conditionally periodic `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_neigh_ncpr` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_neigh_ijcpr(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const ineigh,
    size_t const ndim, size_t const* restrict const nper,
    bool const* const per) {
  size_t jneigh = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = 0; itrial < ktrial; ++itrial)
    if (bmm_neigh_qijcp(pij, ij, itrial, ndim, nper, per)) {
      if (jneigh == ineigh)
        goto br;
      else
        ++jneigh;
    }

  for (size_t itrial = ktrial + 1; itrial < ntrial; ++itrial)
    if (bmm_neigh_qijcp(pij, ij, itrial, ndim, nper, per)) {
      if (jneigh == ineigh)
        goto br;
      else
        ++jneigh;
    }

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    if (per[idim])
      pij[idim] = bmm_size_dec(ij[idim] + pij[idim], 1, nper[idim] - 1);
    else
      pij[idim] = ij[idim] + pij[idim] - 1;
}

/// The call `bmm_neigh_ijcplh(pij, ij, i, d, n, p)`
/// sets the point `pij` to the lower Moore half-neighborhood cell `i`
/// of the point `ij` in a `p`-conditionally periodic `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_neigh_ncplh` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_neigh_ijcplh(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const ineigh,
    size_t const ndim, size_t const* restrict const nper,
    bool const* const per) {
  size_t jneigh = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = 0; itrial < ktrial; ++itrial)
    if (bmm_neigh_qijcp(pij, ij, itrial, ndim, nper, per)) {
      if (jneigh == ineigh)
        goto br;
      else
        ++jneigh;
    }

  if (bmm_neigh_qijcp(pij, ij, ktrial, ndim, nper, per)) {
    if (jneigh == ineigh)
      goto br;
    else
      ++jneigh;
  }

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    if (per[idim])
      pij[idim] = bmm_size_dec(ij[idim] + pij[idim], 1, nper[idim] - 1);
    else
      pij[idim] = ij[idim] + pij[idim] - 1;
}

/// The call `bmm_neigh_ijcplhr(pij, ij, i, d, n, p)`
/// sets the point `pij` to the reduced lower Moore half-neighborhood cell `i`
/// of the point `ij` in a `p`-conditionally periodic `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_neigh_ncplhr` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_neigh_ijcplhr(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const ineigh,
    size_t const ndim, size_t const* restrict const nper,
    bool const* const per) {
  size_t jneigh = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = 0; itrial < ktrial; ++itrial)
    if (bmm_neigh_qijcp(pij, ij, itrial, ndim, nper, per)) {
      if (jneigh == ineigh)
        goto br;
      else
        ++jneigh;
    }

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    if (per[idim])
      pij[idim] = bmm_size_dec(ij[idim] + pij[idim], 1, nper[idim] - 1);
    else
      pij[idim] = ij[idim] + pij[idim] - 1;
}

/// The call `bmm_neigh_ijcpuh(pij, ij, i, d, n, p)`
/// sets the point `pij` to the upper Moore half-neighborhood cell `i`
/// of the point `ij` in a `p`-conditionally periodic `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_neigh_ncpuh` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_neigh_ijcpuh(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const ineigh,
    size_t const ndim, size_t const* restrict const nper,
    bool const* const per) {
  size_t jneigh = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  if (bmm_neigh_qijcp(pij, ij, ktrial, ndim, nper, per)) {
    if (jneigh == ineigh)
      goto br;
    else
      ++jneigh;
  }

  for (size_t itrial = ktrial + 1; itrial < ntrial; ++itrial)
    if (bmm_neigh_qijcp(pij, ij, itrial, ndim, nper, per)) {
      if (jneigh == ineigh)
        goto br;
      else
        ++jneigh;
    }

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    if (per[idim])
      pij[idim] = bmm_size_dec(ij[idim] + pij[idim], 1, nper[idim] - 1);
    else
      pij[idim] = ij[idim] + pij[idim] - 1;
}

/// The call `bmm_neigh_ijcpuhr(pij, ij, i, d, n, p)`
/// sets the point `pij` to the reduced upper Moore half-neighborhood cell `i`
/// of the point `ij` in a `p`-conditionally periodic `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_neigh_ncpuhr` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_neigh_ijcpuhr(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const ineigh,
    size_t const ndim, size_t const* restrict const nper,
    bool const* const per) {
  size_t jneigh = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = ktrial + 1; itrial < ntrial; ++itrial)
    if (bmm_neigh_qijcp(pij, ij, itrial, ndim, nper, per)) {
      if (jneigh == ineigh)
        goto br;
      else
        ++jneigh;
    }

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    if (per[idim])
      pij[idim] = bmm_size_dec(ij[idim] + pij[idim], 1, nper[idim] - 1);
    else
      pij[idim] = ij[idim] + pij[idim] - 1;
}

/// These preprocessor directives allow specifying subsets of neighborhoods.
#define BMM_NEIGH_MASK_NONE (BMM_MASKBITS(1, 0))
#define BMM_NEIGH_MASK_CENTER (BMM_MASKBITS(1, 1))
#define BMM_NEIGH_MASK_LOWER (BMM_MASKBITS(1, 2))
#define BMM_NEIGH_MASK_UPPER (BMM_MASKBITS(1, 3))

// TODO Yes!
/*
__attribute__ ((__nonnull__, __pure__))
inline size_t bmm_neigh_icp(size_t const* restrict const ij,
    size_t const itrial, size_t const ndim,
    size_t const* restrict const nper, bool const* const per,
    int const mask) {
}
*/

/// The call `bmm_neigh_qicp(ij, i, d, n, p)`
/// checks whether the trial index `i` is in the Moore neighborhood
/// of the point `ij` in a `p`-conditionally periodic `n`-wide
/// `d`-dimensional lattice.
/// If there the operation is successful,
/// the corresponding offset index is returned.
/// Otherwise `SIZE_MAX` is returned.
__attribute__ ((__nonnull__, __pure__))
inline size_t bmm_neigh_qicp(size_t const* restrict const ij,
    size_t const itrial, size_t const ndim,
    size_t const* restrict const nper, bool const* const per) {
  size_t icell = 0;

  for (size_t idim = 0; idim < ndim; ++idim) {
    bmm_size_div_t qr = {.quot = itrial, .rem = 0};
    for (size_t jdim = 0; jdim < ndim - idim; ++jdim)
      qr = bmm_size_div(qr.quot, 3);

    if (!per[idim]) {
      size_t const i = ij[idim] + qr.rem;

      if (i <= 0 || i > nper[idim])
        return SIZE_MAX;
    }

    icell *= nper[idim];
    icell += qr.rem;
  }

  return icell;
}

/// The call `bmm_neigh_icpuhr(ij, i, d, n, p)`
/// returns the index of the reduced upper Moore half-neighborhood cell `i`
/// of the point `ij` in a `p`-conditionally periodic `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_neigh_ncpuhr` to find the upper bound of `i`.
__attribute__ ((__nonnull__, __pure__))
inline size_t bmm_neigh_icpuhr(size_t const* restrict const ij,
    size_t const ineigh, size_t const ndim,
    size_t const* restrict const nper, bool const* const per) {
  size_t ioff;
  size_t icell = 0;
  size_t jneigh = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  for (size_t itrial = ktrial + 1; itrial < ntrial; ++itrial) {
    ioff = bmm_neigh_qicp(ij, itrial, ndim, nper, per);

    if (ioff != SIZE_MAX) {
      if (jneigh == ineigh)
        goto br;
      else
        ++jneigh;
    }
  }

br:
  for (size_t idim = 0; idim < ndim; ++idim) {
    bmm_size_div_t qr = {.quot = ioff, .rem = 0};
    for (size_t jdim = 0; jdim < ndim - idim; ++jdim)
      qr = bmm_size_div(qr.quot, nper[ndim - 1 - jdim]);

    if (per[idim])
      qr.rem = bmm_size_dec(ij[idim] + qr.rem, 1, nper[idim] - 1);
    else
      qr.rem = ij[idim] + qr.rem - 1;

    icell *= nper[idim];
    icell += qr.rem;
  }

  return icell;
}

#endif
