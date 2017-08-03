/// Moore neighborhoods that do not overlap themselves.

#ifndef BMM_NEIGH_H
///
/// The naming convention for the procedures in this translation unit
/// obey the following ANTLR 4 grammar.
///
///     grammar Neigh ;
///
///     Query : 'q' ;
///     Number : 'n' ;
///     Iterator : ;
///     Index : 'i' ;
///     IndexVector : 'ij' ;
///     CondPeriodic : 'cp' ;
///     Periodic : 'p' ;
///     Free : ;
///
///     proc : (Query output | Number | Iterator output) bounds input ;
///     input : type ;
///     output : type ;
///     type : Index | IndexVector ;
///     bounds : CondPeriodic | Periodic | Free ;
///
#define BMM_NEIGH_H

#include <alloca.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "common.h"
#include "cpp.h"
#include "ext.h"

/// No neighborhood.
/// This is meant to be used alone and
/// does not need to be included when constructing other masks.
#define BMM_NEIGH_MASK_EMPTY 0

/// The single cell in the center of the neighborhood.
#define BMM_NEIGH_MASK_SINGLE (BMM_MASKBITS(0))

/// The lower half of the neighborhood without the center.
#define BMM_NEIGH_MASK_RLOWERH (BMM_MASKBITS(1))

/// The upper half of the neighborhood without the center.
#define BMM_NEIGH_MASK_RUPPERH (BMM_MASKBITS(2))

/// The lower half of the neighborhood.
#define BMM_NEIGH_MASK_LOWERH (BMM_NEIGH_MASK_SINGLE | BMM_NEIGH_MASK_RLOWERH)

/// The upper half of the neighborhood.
#define BMM_NEIGH_MASK_UPPERH (BMM_NEIGH_MASK_SINGLE | BMM_NEIGH_MASK_RUPPERH)

/// The entire neighborhood.
#define BMM_NEIGH_MASK_FULL (BMM_NEIGH_MASK_LOWERH | BMM_NEIGH_MASK_RUPPERH)

// This implementation is not as optimal as the interface allows,
// but should be good enough for most purposes.

/// The call `bmm_neigh_qijcpij(pijoff, ijcell, iquery, ndim, nper, per)`
/// checks whether the query index `iquery` is in the neighborhood
/// of the cell `ijcell` in a `per`-conditionally periodic `nper`-wide
/// `ndim`-dimensional lattice and
/// sets `pijoff` to the corresponding offset.
__attribute__ ((__nonnull__))
inline bool bmm_neigh_qijcpij(size_t *restrict const pijoff,
    size_t const *restrict const ijcell,
    size_t const iquery, size_t const ndim,
    size_t const *restrict const nper, bool const *const per) {
  type(bmm_hc, size_t)(pijoff, iquery, ndim, 3);

  for (size_t idim = 0; idim < ndim; ++idim) {
    if (per[idim])
      dynamic_assert(nper[idim] >= 5, "Too few neighbor cells");
    else {
      dynamic_assert(nper[idim] >= 3, "Too few neighbor cells");

      size_t const i = ijcell[idim] + pijoff[idim];

      if (i <= 0 || i > nper[idim])
        return false;
    }
  }

  return true;
}

/// The call `bmm_neigh_qijcpi(pijoff, icell, iquery, ndim, nper, per)`
/// checks whether the query index `iquery` is in the neighborhood
/// of the cell `icell` in a `per`-conditionally periodic `nper`-wide
/// `ndim`-dimensional lattice and
/// sets `pijoff` to the corresponding offset.
__attribute__ ((__nonnull__))
inline bool bmm_neigh_qijcpi(size_t *restrict const pijoff,
    size_t const icell, size_t const iquery, size_t const ndim,
    size_t const *restrict const nper, bool const *const per) {
  size_t *const ijcell = alloca(ndim * sizeof *ijcell);

  type(bmm_hcd, size_t)(ijcell, icell, ndim, nper);

  return bmm_neigh_qijcpij(pijoff, ijcell, iquery, ndim, nper, per);
}

/// The call `bmm_neigh_qicpij(ijcell, iquery, ndim, nper, per)`
/// checks whether the query index `iquery` is in the neighborhood
/// of the cell `ijcell` in a `per`-conditionally periodic `nper`-wide
/// `ndim`-dimensional lattice.
/// If there the operation is successful,
/// the corresponding offset is returned.
/// Otherwise `SIZE_MAX` is returned.
__attribute__ ((__nonnull__, __pure__))
inline size_t bmm_neigh_qicpij(size_t const *restrict const ijcell,
    size_t const iquery, size_t const ndim,
    size_t const *restrict const nper, bool const *const per) {
  size_t *const pijoff = alloca(ndim * sizeof *pijoff);

  if (!bmm_neigh_qijcpij(pijoff, ijcell, iquery, ndim, nper, per))
    return SIZE_MAX;

  return type(bmm_unhcd, size_t)(pijoff, ndim, nper);
}

/// The call `bmm_neigh_qicpi(icell, iquery, ndim, nper, per)`
/// checks whether the query index `iquery` is in the neighborhood
/// of the cell `icell` in a `per`-conditionally periodic `nper`-wide
/// `ndim`-dimensional lattice.
/// If there the operation is successful,
/// the corresponding offset is returned.
/// Otherwise `SIZE_MAX` is returned.
__attribute__ ((__nonnull__, __pure__))
inline size_t bmm_neigh_qicpi(size_t const icell,
    size_t const iquery, size_t const ndim,
    size_t const *const nper, bool const *const per) {
  // These temporary variables could be elided with loop fusion,
  // but the result would be quadratic in `ndim`.
  size_t *const pijoff = alloca(ndim * sizeof *pijoff);
  size_t *const ijcell = alloca(ndim * sizeof *ijcell);

  type(bmm_hcd, size_t)(ijcell, icell, ndim, nper);

  if (!bmm_neigh_qijcpij(pijoff, ijcell, iquery, ndim, nper, per))
    return SIZE_MAX;

  return type(bmm_unhcd, size_t)(pijoff, ndim, nper);
}

/// The call `bmm_neigh_ncpij(ijcell, ndim, nper, per, mask)`
/// returns the number of cells in the `mask`-masked neighborhood
/// of the cell `ijcell` in a `per`-conditionally periodic `nper`-wide
/// `ndim`-dimensional lattice.
__attribute__ ((__nonnull__))
inline size_t bmm_neigh_ncpij(size_t const *restrict const ijcell,
    size_t const ndim, size_t const *restrict const nper,
    bool const *const per, int const mask) {
  if (ndim == 0)
    return SIZE_MAX;

  size_t nneigh = 0;

  size_t const nquery = type(bmm_power, size_t)(3, ndim);
  size_t const kquery = nquery / 2;

  if (BMM_MASKANY(mask, BMM_NEIGH_MASK_RLOWERH))
    for (size_t iquery = 0; iquery < kquery; ++iquery)
      if (bmm_neigh_qicpij(ijcell, iquery, ndim, nper, per) != SIZE_MAX)
        ++nneigh;

  if (BMM_MASKANY(mask, BMM_NEIGH_MASK_SINGLE))
    for (size_t iquery = kquery; iquery < kquery + 1; ++iquery)
      if (bmm_neigh_qicpij(ijcell, iquery, ndim, nper, per) != SIZE_MAX)
        ++nneigh;

  if (BMM_MASKANY(mask, BMM_NEIGH_MASK_RUPPERH))
    for (size_t iquery = kquery + 1; iquery < nquery; ++iquery)
      if (bmm_neigh_qicpij(ijcell, iquery, ndim, nper, per) != SIZE_MAX)
        ++nneigh;

  return nneigh;
}

/// The call `bmm_neigh_ncpi(icell, ndim, nper, per, mask)`
/// returns the number of cells in the `mask`-masked neighborhood
/// of the cell `icell` in a `per`-conditionally periodic `nper`-wide
/// `ndim`-dimensional lattice.
__attribute__ ((__nonnull__))
inline size_t bmm_neigh_ncpi(size_t const icell,
    size_t const ndim, size_t const *const nper,
    bool const *const per, int const mask) {
  if (ndim == 0)
    return SIZE_MAX;

  size_t nneigh = 0;

  size_t const nquery = type(bmm_power, size_t)(3, ndim);
  size_t const kquery = nquery / 2;

  if (BMM_MASKANY(mask, BMM_NEIGH_MASK_RLOWERH))
    for (size_t iquery = 0; iquery < kquery; ++iquery)
      if (bmm_neigh_qicpi(icell, iquery, ndim, nper, per) != SIZE_MAX)
        ++nneigh;

  if (BMM_MASKANY(mask, BMM_NEIGH_MASK_SINGLE))
    for (size_t iquery = kquery; iquery < kquery + 1; ++iquery)
      if (bmm_neigh_qicpi(icell, iquery, ndim, nper, per) != SIZE_MAX)
        ++nneigh;

  if (BMM_MASKANY(mask, BMM_NEIGH_MASK_RUPPERH))
    for (size_t iquery = kquery + 1; iquery < nquery; ++iquery)
      if (bmm_neigh_qicpi(icell, iquery, ndim, nper, per) != SIZE_MAX)
        ++nneigh;

  return nneigh;
}

/// The call `bmm_neigh_ijcpij(pijcell, ijcell, ineigh, ndim, nper, per, mask)`
/// sets the cell `pijcell` to the `mask`-masked neighborhood cell `ineigh`
/// of the cell `ijcell` in a `per`-conditionally periodic `nper`-wide
/// `ndim`-dimensional lattice.
/// Use `bmm_neigh_ncpij` to find the upper bound of `ineigh`.
__attribute__ ((__nonnull__))
inline void bmm_neigh_ijcpij(size_t *restrict const pijcell,
    size_t const *restrict const ijcell, size_t const ineigh,
    size_t const ndim, size_t const *restrict const nper,
    bool const *const per, int const mask) {
  size_t const nquery = type(bmm_power, size_t)(3, ndim);
  size_t const kquery = nquery / 2;

  size_t jneigh = 0;

  if (BMM_MASKANY(mask, BMM_NEIGH_MASK_RLOWERH))
    for (size_t iquery = 0; iquery < kquery; ++iquery)
      if (bmm_neigh_qijcpij(pijcell, ijcell, iquery, ndim, nper, per)) {
        if (jneigh == ineigh)
          goto br;
        else
          ++jneigh;
      }

  if (BMM_MASKANY(mask, BMM_NEIGH_MASK_SINGLE))
    for (size_t iquery = kquery; iquery < kquery + 1; ++iquery)
      if (bmm_neigh_qijcpij(pijcell, ijcell, iquery, ndim, nper, per)) {
        if (jneigh == ineigh)
          goto br;
        else
          ++jneigh;
      }

  if (BMM_MASKANY(mask, BMM_NEIGH_MASK_RUPPERH))
    for (size_t iquery = kquery + 1; iquery < nquery; ++iquery)
      if (bmm_neigh_qijcpij(pijcell, ijcell, iquery, ndim, nper, per)) {
        if (jneigh == ineigh)
          goto br;
        else
          ++jneigh;
      }

  return;

br:
  for (size_t idim = 0; idim < ndim; ++idim)
    if (per[idim])
      pijcell[idim] = type(bmm_dec, size_t)(ijcell[idim] + pijcell[idim],
          1, nper[idim] - 1);
    else
      pijcell[idim] = ijcell[idim] + pijcell[idim] - 1;
}

/// The call `bmm_neigh_ijcpi(pijcell, icell, ineigh, ndim, nper, per, mask)`
/// sets the cell `pijcell` to the `mask`-masked neighborhood cell `ineigh`
/// of the cell `icell` in a `per`-conditionally periodic `nper`-wide
/// `ndim`-dimensional lattice.
/// Use `bmm_neigh_ncpij` to find the upper bound of `ineigh`.
__attribute__ ((__nonnull__))
inline void bmm_neigh_ijcpi(size_t *restrict const pijcell,
    size_t const icell, size_t const ineigh,
    size_t const ndim, size_t const *restrict const nper,
    bool const *const per, int const mask) {
  size_t *const ijcell = alloca(ndim * sizeof *ijcell);

  type(bmm_hcd, size_t)(ijcell, icell, ndim, nper);

  bmm_neigh_ijcpij(pijcell, ijcell, ineigh, ndim, nper, per, mask);
}

/// The call `bmm_neigh_icpij(ijcell, ineigh, ndim, nper, per, mask)`
/// returns the index of the `mask`-masked neighborhood cell `ineigh`
/// of the cell `ijcell` in a `per`-conditionally periodic `nper`-wide
/// `ndim`-dimensional lattice.
/// Use `bmm_neigh_ncpij` to find the upper bound of `ineigh`.
__attribute__ ((__nonnull__, __pure__))
inline size_t bmm_neigh_icpij(size_t const *restrict const ijcell,
    size_t const ineigh, size_t const ndim,
    size_t const *restrict const nper, bool const *const per,
    int const mask) {
  size_t *const pijcell = alloca(ndim * sizeof *pijcell);

  bmm_neigh_ijcpij(pijcell, ijcell, ineigh, ndim, nper, per, mask);

  return type(bmm_unhcd, size_t)(pijcell, ndim, nper);
}

/// The call `bmm_neigh_icpi(icell, ineigh, ndim, nper, per, mask)`
/// returns the index of the `mask`-masked neighborhood cell `ineigh`
/// of the cell `icell` in a `per`-conditionally periodic `nper`-wide
/// `ndim`-dimensional lattice.
/// Use `bmm_neigh_ncpi` to find the upper bound of `icell`.
__attribute__ ((__nonnull__, __pure__))
inline size_t bmm_neigh_icpi(size_t const icell,
    size_t const ineigh, size_t const ndim,
    size_t const *restrict const nper, bool const *const per,
    int const mask) {
  size_t *const pijcell = alloca(ndim * sizeof *pijcell);
  size_t *const ijcell = alloca(ndim * sizeof *ijcell);

  type(bmm_hcd, size_t)(ijcell, icell, ndim, nper);

  bmm_neigh_ijcpij(pijcell, ijcell, ineigh, ndim, nper, per, mask);

  return type(bmm_unhcd, size_t)(pijcell, ndim, nper);
}

#endif
