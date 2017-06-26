#ifndef BMM_NEIGH_H
/// Moore neighborhoods.
#define BMM_NEIGH_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "cpp.h"
#include "ext.h"
#include "size.h"

// This implementation is not as optimal as the interface allows,
// but should be good enough for most purposes.

/// No neighborhood.
/// This does not need to be included when constructing other masks.
#define BMM_NEIGH_MASK_EMPTY 0

/// The single cell in the center of the neighborhood.
#define BMM_NEIGH_MASK_SINGLE (BMM_MASKBITS(1, 0))

/// The lower half of the neighborhood without the center.
#define BMM_NEIGH_MASK_RLOWERH (BMM_MASKBITS(1, 1))

/// The upper half of the neighborhood without the center.
#define BMM_NEIGH_MASK_RUPPERH (BMM_MASKBITS(1, 2))

/// The lower half of the neighborhood.
#define BMM_NEIGH_MASK_LOWERH (BMM_NEIGH_MASK_SINGLE | BMM_NEIGH_MASK_RLOWERH)

/// The upper half of the neighborhood.
#define BMM_NEIGH_MASK_UPPERH (BMM_NEIGH_MASK_SINGLE | BMM_NEIGH_MASK_RUPPERH)

/// The entire neighborhood.
#define BMM_NEIGH_MASK_FULL (BMM_NEIGH_MASK_LOWERH | BMM_NEIGH_MASK_RUPPERH)

/// The call `bmm_neigh_qijcp(pij, ij, i, d, n, p)`
/// checks whether the trial index `i` is in the neighborhood
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

/// The call `bmm_neigh_qicp(ij, i, d, n, p)`
/// checks whether the trial index `i` is in the neighborhood
/// of the point `ij` in a `p`-conditionally periodic `n`-wide
/// `d`-dimensional lattice.
/// If there the operation is successful,
/// the corresponding offset is returned.
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

// TODO Finish this test.
/// The call `bmm_neigh_qicplin(j, i, d, n, p)`
/// checks whether the trial index `i` is in the neighborhood
/// of the point `j` in a `p`-conditionally periodic `n`-wide
/// `d`-dimensional lattice.
/// If there the operation is successful,
/// the corresponding offset is returned.
/// Otherwise `SIZE_MAX` is returned.
__attribute__ ((__nonnull__, __pure__))
inline size_t bmm_neigh_qicplin(size_t const i,
    size_t const itrial, size_t const ndim,
    size_t const* const nper, bool const* const per) {
  size_t icell = 0;

  for (size_t idim = 0; idim < ndim; ++idim) {
    bmm_size_div_t qrlin = {.quot = i, .rem = 0};
    bmm_size_div_t qr = {.quot = itrial, .rem = 0};
    for (size_t jdim = 0; jdim < ndim - idim; ++jdim) {
      qrlin = bmm_size_div(qrlin.quot, nper[ndim - 1 - jdim]);
      qr = bmm_size_div(qr.quot, 3);
    }

    if (!per[idim]) {
      size_t const i = qrlin.rem + qr.rem;

      if (i <= 0 || i > nper[idim])
        return SIZE_MAX;
    }

    icell *= nper[idim];
    icell += qr.rem;
  }

  return icell;
}

/// The call `bmm_neigh_ncp(ij, d, n, p, m)`
/// returns the number of cells in the `m`-masked neighborhood
/// of the point `ij` in a `p`-conditionally periodic `n`-wide
/// `d`-dimensional lattice.
__attribute__ ((__nonnull__))
inline size_t bmm_neigh_ncp(size_t const* restrict const ij,
    size_t const ndim, size_t const* restrict const nper,
    bool const* const per, int const mask) {
  if (ndim == 0)
    return SIZE_MAX;

  size_t nneigh = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  if (BMM_MASKANY(mask, BMM_NEIGH_MASK_RLOWERH))
    for (size_t itrial = 0; itrial < ktrial; ++itrial)
      if (bmm_neigh_qicp(ij, itrial, ndim, nper, per) != SIZE_MAX)
        ++nneigh;

  if (BMM_MASKANY(mask, BMM_NEIGH_MASK_SINGLE))
    for (size_t itrial = ktrial; itrial < ktrial + 1; ++itrial)
      if (bmm_neigh_qicp(ij, itrial, ndim, nper, per) != SIZE_MAX)
        ++nneigh;

  if (BMM_MASKANY(mask, BMM_NEIGH_MASK_RUPPERH))
    for (size_t itrial = ktrial + 1; itrial < ntrial; ++itrial)
      if (bmm_neigh_qicp(ij, itrial, ndim, nper, per) != SIZE_MAX)
        ++nneigh;

  return nneigh;
}

/// The call `bmm_neigh_ijcp(pij, ij, i, d, n, p, m)`
/// sets the point `pij` to the `m`-masked neighborhood cell `i`
/// of the point `ij` in a `p`-conditionally periodic `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_neigh_ncpuhr` to find the upper bound of `i`.
__attribute__ ((__nonnull__))
inline void bmm_neigh_ijcp(size_t* restrict const pij,
    size_t const* restrict const ij, size_t const ineigh,
    size_t const ndim, size_t const* restrict const nper,
    bool const* const per, int const mask) {
  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  size_t jneigh = 0;

  if (BMM_MASKANY(mask, BMM_NEIGH_MASK_RLOWERH))
    for (size_t itrial = 0; itrial < ktrial; ++itrial)
      if (bmm_neigh_qijcp(pij, ij, itrial, ndim, nper, per)) {
        if (jneigh == ineigh)
          goto br;
        else
          ++jneigh;
      }

  if (BMM_MASKANY(mask, BMM_NEIGH_MASK_SINGLE))
    for (size_t itrial = ktrial; itrial < ktrial + 1; ++itrial)
      if (bmm_neigh_qijcp(pij, ij, itrial, ndim, nper, per)) {
        if (jneigh == ineigh)
          goto br;
        else
          ++jneigh;
      }

  if (BMM_MASKANY(mask, BMM_NEIGH_MASK_RUPPERH))
    for (size_t itrial = ktrial + 1; itrial < ntrial; ++itrial)
      if (bmm_neigh_qijcp(pij, ij, itrial, ndim, nper, per)) {
        if (jneigh == ineigh)
          goto br;
        else
          ++jneigh;
      }

  return;

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    if (per[idim])
      pij[idim] = bmm_size_dec(ij[idim] + pij[idim], 1, nper[idim] - 1);
    else
      pij[idim] = ij[idim] + pij[idim] - 1;
}

/// The call `bmm_neigh_icp(ij, i, d, n, p, m)`
/// returns the index of the `m`-masked neighborhood cell `i`
/// of the point `ij` in a `p`-conditionally periodic `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_neigh_ncp` to find the upper bound of `i`.
__attribute__ ((__nonnull__, __pure__))
inline size_t bmm_neigh_icp(size_t const* restrict const ij,
    size_t const ineigh, size_t const ndim,
    size_t const* restrict const nper, bool const* const per,
    int const mask) {
  size_t ioff;
  size_t icell = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  size_t jneigh = 0;

  if (BMM_MASKANY(mask, BMM_NEIGH_MASK_RLOWERH))
    for (size_t itrial = 0; itrial < ktrial; ++itrial) {
      ioff = bmm_neigh_qicp(ij, itrial, ndim, nper, per);

      if (ioff != SIZE_MAX) {
        if (jneigh == ineigh)
          goto br;
        else
          ++jneigh;
      }
    }

  if (BMM_MASKANY(mask, BMM_NEIGH_MASK_SINGLE))
    for (size_t itrial = ktrial; itrial < ktrial + 1; ++itrial) {
      ioff = bmm_neigh_qicp(ij, itrial, ndim, nper, per);

      if (ioff != SIZE_MAX) {
        if (jneigh == ineigh)
          goto br;
        else
          ++jneigh;
      }
    }

  if (BMM_MASKANY(mask, BMM_NEIGH_MASK_RUPPERH))
    for (size_t itrial = ktrial + 1; itrial < ntrial; ++itrial) {
      ioff = bmm_neigh_qicp(ij, itrial, ndim, nper, per);

      if (ioff != SIZE_MAX) {
        if (jneigh == ineigh)
          goto br;
        else
          ++jneigh;
      }
    }

  return SIZE_MAX;

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

// TODO Finish this test.
/// The call `bmm_neigh_icplin(j, i, d, n, p, m)`
/// returns the index of the `m`-masked neighborhood cell `i`
/// of the point `j` in a `p`-conditionally periodic `n`-wide
/// `d`-dimensional lattice.
/// Use `bmm_neigh_ncp` to find the upper bound of `i`.
__attribute__ ((__nonnull__, __pure__))
inline size_t bmm_neigh_icplin(size_t const i,
    size_t const ineigh, size_t const ndim,
    size_t const* restrict const nper, bool const* const per,
    int const mask) {
  size_t ioff;
  size_t icell = 0;

  size_t const ntrial = bmm_size_pow(3, ndim);
  size_t const ktrial = ntrial / 2;

  size_t jneigh = 0;

  if (BMM_MASKANY(mask, BMM_NEIGH_MASK_RLOWERH))
    for (size_t itrial = 0; itrial < ktrial; ++itrial) {
      ioff = bmm_neigh_qicplin(i, itrial, ndim, nper, per);

      if (ioff != SIZE_MAX) {
        if (jneigh == ineigh)
          goto br;
        else
          ++jneigh;
      }
    }

  if (BMM_MASKANY(mask, BMM_NEIGH_MASK_SINGLE))
    for (size_t itrial = ktrial; itrial < ktrial + 1; ++itrial) {
      ioff = bmm_neigh_qicplin(i, itrial, ndim, nper, per);

      if (ioff != SIZE_MAX) {
        if (jneigh == ineigh)
          goto br;
        else
          ++jneigh;
      }
    }

  if (BMM_MASKANY(mask, BMM_NEIGH_MASK_RUPPERH))
    for (size_t itrial = ktrial + 1; itrial < ntrial; ++itrial) {
      ioff = bmm_neigh_qicplin(i, itrial, ndim, nper, per);

      if (ioff != SIZE_MAX) {
        if (jneigh == ineigh)
          goto br;
        else
          ++jneigh;
      }
    }

  return SIZE_MAX;

br:
  for (size_t idim = 0; idim < ndim; ++idim) {
    bmm_size_div_t qrlin = {.quot = i, .rem = 0};
    bmm_size_div_t qr = {.quot = ioff, .rem = 0};
    for (size_t jdim = 0; jdim < ndim - idim; ++jdim) {
      qrlin = bmm_size_div(qrlin.quot, nper[ndim - 1 - jdim]);
      qr = bmm_size_div(qr.quot, nper[ndim - 1 - jdim]);
    }

    if (per[idim])
      qr.rem = bmm_size_dec(qrlin.rem + qr.rem, 1, nper[idim] - 1);
    else
      qr.rem = qrlin.rem + qr.rem - 1;

    icell *= nper[idim];
    icell += qr.rem;
  }

  return icell;
}

#endif
