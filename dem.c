#include <gsl/gsl_rng.h>
#include <math.h>
#include <signal.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifdef _GNU_SOURCE
#include <fenv.h>
#endif

#include "common.h"
#include "conf.h"
#include "cpp.h"
#include "dem.h"
#include "fp.h"
#include "geom.h"
#include "geom2d.h"
#include "io.h"
#include "ival.h"
#include "kernel.h"
#include "msg.h"
#include "neigh.h"
#include "random.h"
#include "sig.h"
#include "tle.h"

__attribute__ ((__nonnull__))
double bmm_dem_est_vdc(double *const pv, struct bmm_dem const *const dem) {
  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    pv[idim] = 0.0;

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    if (dem->part.role[ipart] == BMM_DEM_ROLE_DRIVEN)
      for (size_t idim = 0; idim < BMM_NDIM; ++idim)
        pv[idim] += dem->part.v[ipart][idim];

  size_t ndrive = 0;
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    if (dem->part.role[ipart] == BMM_DEM_ROLE_DRIVEN)
      ++ndrive;

  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    pv[idim] /= (double) ndrive;
}

size_t bmm_dem_script_addstage(struct bmm_dem_opts *const opts) {
  size_t const istage = opts->script.n;

  if (istage >= BMM_MSTAGE)
    return SIZE_MAX;

  ++opts->script.n;

  opts->script.mode[istage] = BMM_DEM_MODE_IDLE;
  opts->script.tspan[istage] = 0.0;
  opts->script.dt[istage] = 0.0;

  return istage;
}

bool bmm_dem_script_ongoing(struct bmm_dem const *const dem) {
  return dem->script.i < dem->opts.script.n;
}

bool bmm_dem_script_trans(struct bmm_dem *const dem) {
  double const toff = dem->time.t - dem->script.tprev -
    dem->opts.script.tspan[dem->script.i];

  if (toff >= 0.0) {
    dem->script.tprev = dem->time.t;
    dem->script.ttrans[dem->script.i] = dem->time.t;
    dem->script.toff[dem->script.i] = toff;
    ++dem->script.i;

    return bmm_dem_script_ongoing(dem);
  }

  return true;
}

/// The call `bmm_dem_cache_j(dem, ipart)`
/// caches the moment of inertia of the particle `ipart`
/// in the simulation `dem`.
__attribute__ ((__nonnull__))
static void bmm_dem_cache_j(struct bmm_dem *const dem,
    size_t const ipart) {
  dem->cache.j[ipart] = dem->part.jred[ipart] *
    dem->part.m[ipart] * $(bmm_power, double)(dem->part.r[ipart], 2);
}

/// The call `bmm_dem_cache_x(dem, ipart)`
/// caches the position of the particle `ipart`
/// in the simulation `dem`.
__attribute__ ((__nonnull__))
static void bmm_dem_cache_x(struct bmm_dem *const dem,
    size_t const ipart) {
  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    dem->cache.x[ipart][idim] = dem->part.x[ipart][idim];
}

/// The call `bmm_dem_cache_ijcell(dem, ipart)`
/// caches the neighbor cell index vector of the particle `ipart`
/// in the simulation `dem`.
/// The positions need to be cached first
/// by calling `bmm_dem_cache_x`.
__attribute__ ((__nonnull__))
static void bmm_dem_cache_ijcell(struct bmm_dem *const dem,
    size_t const ipart) {
  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    dem->cache.ijcell[ipart][idim] = bmm_fp_iclerp(dem->cache.x[ipart][idim],
        0.0, dem->opts.box.x[idim], 1, dem->opts.cache.ncell[idim] - 1);
}

/// The call `bmm_dem_cache_icell(dem, ipart)`
/// caches the neighbor cell index of the particle `ipart`
/// in the simulation `dem`.
/// The neighbor cell index vectors need to be cached first
/// by calling `bmm_dem_cache_ijcell`.
__attribute__ ((__nonnull__))
static void bmm_dem_cache_icell(struct bmm_dem *const dem,
    size_t const ipart) {
  dem->cache.icell[ipart] = $(bmm_unhcd, size_t)(dem->cache.ijcell[ipart],
      BMM_NDIM, dem->opts.cache.ncell);
}

/// The call `bmm_dem_cache_clrparts(dem)`
/// clears the neighbor cell index mapping cache
/// in the simulation `dem`.
__attribute__ ((__nonnull__))
static void bmm_dem_cache_clrparts(struct bmm_dem *const dem) {
  for (size_t icell = 0; icell < nmembof(dem->cache.part); ++icell)
    dem->cache.part[icell].n = 0;
}

/// The call `bmm_dem_cache_addpart(dem, ipart)`
/// tries to cache the neighbor cell index mapping of the particle `ipart`
/// in the simulation `dem`.
/// If there is enough capacity and the operation is successful,
/// `true` is returned.
/// Otherwise `false` is returned and
/// the previous neighbor relations are restored.
/// The neighbor cell indices need to be cached first
/// by calling `bmm_dem_cache_icell`.
__attribute__ ((__nonnull__))
static bool bmm_dem_cache_addpart(struct bmm_dem *const dem,
    size_t const ipart) {
  size_t const icell = dem->cache.icell[ipart];

  if (dem->cache.part[icell].n >= nmembof(dem->cache.part[icell].i)) {
    BMM_TLE_EXTS(BMM_TLE_NUM_UNKNOWN, "Neighbor cell capacity exceeded");

    return false;
  }

  dem->cache.part[icell].i[dem->cache.part[icell].n] = ipart;
  ++dem->cache.part[icell].n;

  return true;
}

/// The call `bmm_dem_cache_eligible(dem, ipart, jpart)`
/// checks whether the particles `ipart` and `jpart` are eligible neighbors
/// in the simulation `dem`.
/// The neighbor cell indices need to be cached first
/// by calling `bmm_dem_cache_icell`.
/// Note that this function is neither symmetric nor reflexive
/// with respect to particle indices.
__attribute__ ((__nonnull__, __pure__))
static bool bmm_dem_cache_eligible(struct bmm_dem const *const dem,
    size_t const ipart, size_t const jpart) {
  if (dem->cache.icell[ipart] == dem->cache.icell[jpart] && jpart <= ipart)
    return false;

  if (bmm_geom2d_cpdist2(dem->cache.x[ipart], dem->cache.x[jpart],
        dem->opts.box.x, dem->opts.box.per) >
      $(bmm_power, double)(dem->opts.cache.dcutoff, 2))
    return false;

  return true;
}

/// The call `bmm_dem_cache_clrneighs(dem)`
/// clears the neighbor cache
/// in the simulation `dem`.
__attribute__ ((__nonnull__))
static void bmm_dem_cache_clrneighs(struct bmm_dem *const dem) {
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    dem->cache.neigh[ipart].n = 0;
}

// The procedures `bmm_dem_cache_addfrom` and `bmm_dem_cache_addto`
// are identical except for those parts that are marked
// covariant or contravariant with respect to particle indices.

/// The call `bmm_dem_cache_addfrom(dem, ipart, mask)`
/// tries to add all the eligible particles
/// inside the `mask`-masked neighborhood
/// of the particle `ipart` to its neighbors
/// in the simulation `dem`.
/// If there is enough capacity and the operation is successful,
/// `true` is returned.
/// Otherwise `false` is returned and
/// the cache is left in an undefined state.
/// The neighbor cell index mappings need to be cached first
/// by calling `bmm_dem_cache_addpart`.
__attribute__ ((__nonnull__))
static bool bmm_dem_cache_addfrom(struct bmm_dem *const dem,
    size_t const ipart, int const mask) {
  size_t const nneigh = bmm_neigh_ncpij(dem->cache.ijcell[ipart],
      BMM_NDIM, dem->opts.cache.ncell, dem->opts.box.per, mask);

  for (size_t ineigh = 0; ineigh < nneigh; ++ineigh) {
    size_t const icell = bmm_neigh_icpij(dem->cache.ijcell[ipart], ineigh,
        BMM_NDIM, dem->opts.cache.ncell, dem->opts.box.per, mask);

    for (size_t igroup = 0; igroup < dem->cache.part[icell].n; ++igroup) {
      size_t const jpart = dem->cache.part[icell].i[igroup];

      // This block is covariant.
      if (bmm_dem_cache_eligible(dem, ipart, jpart)) {
        if (dem->cache.neigh[ipart].n >= nmembof(dem->cache.neigh[ipart].i)) {
          BMM_TLE_EXTS(BMM_TLE_NUM_UNKNOWN, "Neighbor capacity exceeded");

          return false;
        }

        dem->cache.neigh[ipart].i[dem->cache.neigh[ipart].n] = jpart;
        ++dem->cache.neigh[ipart].n;
      }
    }
  }

  return true;
}

/// The call `bmm_dem_cache_addto(dem, ipart)`
/// tries to add the particle `ipart` to the neighbors
/// of all the eligible particles
/// inside its `mask`-masked neighborhood
/// in the simulation `dem`.
/// If there is enough capacity and the operation is successful,
/// `true` is returned.
/// Otherwise `false` is returned and
/// the cache is left in an undefined state.
/// The neighbor cell index mappings need to be cached first
/// by calling `bmm_dem_cache_addpart`.
__attribute__ ((__nonnull__))
static bool bmm_dem_cache_addto(struct bmm_dem *const dem,
    size_t const ipart, int const mask) {
  size_t const nneigh = bmm_neigh_ncpij(dem->cache.ijcell[ipart],
      BMM_NDIM, dem->opts.cache.ncell, dem->opts.box.per, mask);

  for (size_t ineigh = 0; ineigh < nneigh; ++ineigh) {
    size_t const icell = bmm_neigh_icpij(dem->cache.ijcell[ipart], ineigh,
        BMM_NDIM, dem->opts.cache.ncell, dem->opts.box.per, mask);

    for (size_t igroup = 0; igroup < dem->cache.part[icell].n; ++igroup) {
      size_t const jpart = dem->cache.part[icell].i[igroup];

      // This block is contravariant.
      if (bmm_dem_cache_eligible(dem, jpart, ipart)) {
        if (dem->cache.neigh[jpart].n >= nmembof(dem->cache.neigh[jpart].i)) {
          BMM_TLE_EXTS(BMM_TLE_NUM_UNKNOWN, "Neighbor capacity exceeded");

          return false;
        }

        dem->cache.neigh[jpart].i[dem->cache.neigh[jpart].n] = ipart;
        ++dem->cache.neigh[jpart].n;
      }
    }
  }

  return true;
}

bool bmm_dem_cache_build(struct bmm_dem *const dem) {
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    bmm_dem_cache_j(dem, ipart);

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    bmm_dem_cache_x(dem, ipart);

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    bmm_dem_cache_ijcell(dem, ipart);

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    bmm_dem_cache_icell(dem, ipart);

  bmm_dem_cache_clrparts(dem);
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    if (!bmm_dem_cache_addpart(dem, ipart))
      return false;

  bmm_dem_cache_clrneighs(dem);
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    if (!bmm_dem_cache_addfrom(dem, ipart, BMM_NEIGH_MASK_UPPERH))
      return false;

  dem->cache.stale = false;

  return true;
}

size_t bmm_dem_search_cont(struct bmm_dem const *const dem,
    enum bmm_dem_ct const ict, size_t const ipart, size_t const jpart) {
  for (size_t icont = 0; icont < dem->pair[ict].cont.src[ipart].n; ++icont) {
    size_t const kpart = dem->pair[ict].cont.src[ipart].itgt[icont];

    if (kpart == jpart)
      return icont;
  }

  return SIZE_MAX;
}

size_t bmm_dem_addcont_unsafe(struct bmm_dem *const dem,
    enum bmm_dem_ct const ict, size_t const ipart, size_t const jpart) {
  size_t const icont = dem->pair[ict].cont.src[ipart].n;

  if (icont >= BMM_MCONTACT)
    return SIZE_MAX;

  ++dem->pair[ict].cont.src[ipart].n;

  double xdiffij[BMM_NDIM];
  bmm_geom2d_cpdiff(xdiffij, dem->part.x[jpart], dem->part.x[ipart],
      dem->opts.box.x, dem->opts.box.per);

  double const d = bmm_geom2d_norm(xdiffij);
  double const r = dem->part.r[ipart] + dem->part.r[jpart];

  // TODO Excessive wrapping!
  double const gammaij = $(bmm_swrap, double)(bmm_geom2d_dir(xdiffij), M_2PI);
  double const gammaji = $(bmm_swrap, double)(gammaij + M_PI, M_2PI);
  double const phii = $(bmm_swrap, double)(dem->part.phi[ipart], M_2PI);
  double const phij = $(bmm_swrap, double)(dem->part.phi[jpart], M_2PI);
  double const chii = $(bmm_swrap, double)(gammaij - phii, M_2PI);
  double const chij = $(bmm_swrap, double)(gammaji - phij, M_2PI);

  dem->pair[ict].cont.src[ipart].chirest[icont][0] = chii;

  dem->pair[ict].cont.src[ipart].chirest[icont][1] = chij;

  double a = 0.3;
  double v[2];
  v[0] = 1.0 - a;
  v[1] = 1.0 + a;
  dem->pair[ict].cont.src[ipart].strength[icont] = bmm_random_get(dem->rng, v);

  double const drest = d * dem->opts.cont.cshcont;

  dem->pair[ict].cont.src[ipart].drest[icont] = drest;
  dem->pair[ict].cont.src[ipart].itgt[icont] = jpart;

  dem->pair[ict].cont.src[ipart].epsilon[icont] = 0.0;

  // TODO Really?
  // dem->cache.stale = true;

  return icont;
}

size_t bmm_dem_addcont(struct bmm_dem *const dem,
    enum bmm_dem_ct const ict, size_t const ipart, size_t const jpart) {
  // TODO Does this ever fail?
  // dynamic_assert(ipart < jpart, "Yes");
  if (ipart < jpart)
    return bmm_dem_addcont_unsafe(dem, ict, ipart, jpart);
  else
    return bmm_dem_addcont_unsafe(dem, ict, jpart, ipart);
}

/// The call `bmm_dem_copycont(dem, ipart, icont, jcont)`
/// reassigns the contact `jcont` to `icont` for the particle `ipart`.
__attribute__ ((__nonnull__))
static void bmm_dem_copycont(struct bmm_dem *const dem,
    enum bmm_dem_ct const ict, size_t const ipart, size_t const icont, size_t const jcont) {
  dem->pair[ict].cont.src[ipart].itgt[icont] = dem->pair[ict].cont.src[ipart].itgt[jcont];

  dem->pair[ict].cont.src[ipart].drest[icont] = dem->pair[ict].cont.src[ipart].drest[jcont];

  for (size_t iend = 0; iend < BMM_NEND; ++iend)
    dem->pair[ict].cont.src[ipart].chirest[icont][iend] =
      dem->pair[ict].cont.src[ipart].chirest[jcont][iend];

  dem->pair[ict].cont.src[ipart].rlim[icont] = dem->pair[ict].cont.src[ipart].rlim[jcont];

  dem->pair[ict].cont.src[ipart].philim[icont] = dem->pair[ict].cont.src[ipart].philim[jcont];

  dem->pair[ict].cont.src[ipart].epsilon[icont] = dem->pair[ict].cont.src[ipart].epsilon[jcont];
}

void bmm_dem_remcont_unsafe(struct bmm_dem *const dem,
    enum bmm_dem_ct const ict, size_t const ipart, size_t const icont) {
  --dem->pair[ict].cont.src[ipart].n;

  size_t const jcont = dem->pair[ict].cont.src[ipart].n;

  if (jcont != icont)
    bmm_dem_copycont(dem, ict, ipart, icont, jcont);

  // TODO Really?
  // dem->cache.stale = true;
}

void bmm_dem_remcont(struct bmm_dem *const dem,
    enum bmm_dem_ct const ict, size_t const ipart, size_t const jpart, size_t const icont) {
  // TODO Does this ever fail?
  // dynamic_assert(ipart < jpart, "Yes");
  if (ipart < jpart)
    bmm_dem_remcont_unsafe(dem, ict, ipart, icont);
  else
    bmm_dem_remcont_unsafe(dem, ict, jpart, icont);
}

bool bmm_dem_yield_pair(struct bmm_dem *const dem,
    size_t const ipart, size_t const jpart, size_t const icont) {
  double xdiffji[BMM_NDIM];
  bmm_geom2d_cpdiff(xdiffji, dem->part.x[ipart], dem->part.x[jpart],
      dem->opts.box.x, dem->opts.box.per);

  double const d2 = bmm_geom2d_norm2(xdiffji);
  double const d = sqrt(d2);

  /*
  // TODO Do not use `rlim` but rather `f` and proper stress criteria.
  if (d2 > $(bmm_power, double)(dem->pair[BMM_DEM_CT_STRONG].cont.src[ipart].rlim[icont], 2)) {
    bmm_dem_remcont(dem, BMM_DEM_CT_STRONG, ipart, jpart, icont);

    return true;
  }
  */

  double xnormji[BMM_NDIM];
  bmm_geom2d_scale(xnormji, xdiffji, 1.0 / d);

  double xtangji[BMM_NDIM];
  bmm_geom2d_rperp(xtangji, xnormji);

  double fdiffij[BMM_NDIM];
  bmm_geom2d_diff(fdiffij, dem->part.f[jpart], dem->part.f[ipart]);

  double const fnormij = bmm_geom2d_dot(fdiffij, xnormji);
  double const ftangij = $(bmm_abs, double)(bmm_geom2d_dot(fdiffij, xtangji));

  // TODO Is this a valid assumption for the yield point cross section?
  double const a = $(bmm_hmean2, double)(dem->part.r[ipart], dem->part.r[jpart]);

  // TODO Does this make sense here?
  double const strength = dem->pair[BMM_DEM_CT_STRONG].cont.src[ipart].strength[icont];

  double const sigmaij = fnormij / a / strength;
  double const tauij = ftangij / a / strength;

  switch (dem->yield.tag) {
    case BMM_DEM_YIELD_RANKINE:
      if (sigmaij >= dem->yield.params.rankine.sigmacrit) {
        bmm_dem_remcont(dem, BMM_DEM_CT_STRONG, ipart, jpart, icont);

        return true;
      }

      break;
    case BMM_DEM_YIELD_TRESCA:
      if (sigmaij >= dem->yield.params.tresca.taucrit) {
        bmm_dem_remcont(dem, BMM_DEM_CT_STRONG, ipart, jpart, icont);

        return true;
      }

      break;
    case BMM_DEM_YIELD_ZE:
      if ($(bmm_power, double)(sigmaij / dem->yield.params.ze.sigmacrit, 2) +
          $(bmm_power, double)(tauij / dem->yield.params.ze.taucrit, 2) >= 1.0) {
        bmm_dem_remcont(dem, BMM_DEM_CT_STRONG, ipart, jpart, icont);

        return true;
      }

      break;
  }

  return false;
}

void bmm_dem_analyze_pair(struct bmm_dem *const dem,
    size_t const ipart, size_t const jpart) {
  // TODO Settle these order problems.
  if (ipart >= jpart) {
    bmm_dem_analyze_pair(dem, jpart, ipart);

    return;
  }

  size_t const icont = bmm_dem_search_cont(dem, BMM_DEM_CT_STRONG, ipart, jpart);
  if (icont != SIZE_MAX)
    (void) bmm_dem_yield_pair(dem, ipart, jpart, icont);
  else {
    double xdiffji[BMM_NDIM];
    bmm_geom2d_cpdiff(xdiffji, dem->part.x[ipart], dem->part.x[jpart],
        dem->opts.box.x, dem->opts.box.per);

    double const d2 = bmm_geom2d_norm2(xdiffji);
    double const r2 = $(bmm_power, double)(dem->part.r[ipart] + dem->part.r[jpart], 2);
    bool const overlap = d2 <= r2;

    size_t const jcont = bmm_dem_search_cont(dem, BMM_DEM_CT_WEAK, ipart, jpart);
    if (jcont != SIZE_MAX) {
      if (!overlap)
        bmm_dem_remcont(dem, BMM_DEM_CT_WEAK, ipart, jpart, jcont);
    } else {
      if (overlap)
        bmm_dem_addcont(dem, BMM_DEM_CT_WEAK, ipart, jpart);
    }
  }
}

void bmm_dem_analyze(struct bmm_dem *const dem) {
  switch (dem->cache.tag) {
    case BMM_DEM_CACHE_NONE:
      for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
        for (size_t jpart = ipart + 1; jpart < dem->part.n; ++jpart)
          bmm_dem_analyze_pair(dem, ipart, jpart);

      break;
    case BMM_DEM_CACHE_NEIGH:
      for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
        for (size_t ineigh = 0; ineigh < dem->cache.neigh[ipart].n; ++ineigh) {
          size_t const jpart = dem->cache.neigh[ipart].i[ineigh];

          bmm_dem_analyze_pair(dem, ipart, jpart);
        }

      break;
  }
}

size_t bmm_dem_addpart(struct bmm_dem *const dem) {
  size_t const ipart = dem->part.n;

  if (ipart >= BMM_MPART)
    return SIZE_MAX;

  ++dem->part.n;
  dem->part.l[ipart] = dem->part.lnew;
  ++dem->part.lnew;

  dem->part.role[ipart] = BMM_DEM_ROLE_FREE;
  dem->part.r[ipart] = 1.0;
  dem->part.m[ipart] = 1.0;
  dem->part.jred[ipart] = bmm_geom_ballprmoi(BMM_NDIM);

  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    dem->part.x[ipart][idim] = 0.0;

  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    dem->part.v[ipart][idim] = 0.0;

  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    dem->part.a[ipart][idim] = 0.0;

  dem->part.phi[ipart] = 0.0;
  dem->part.omega[ipart] = 0.0;
  dem->part.alpha[ipart] = 0.0;

  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    dem->part.f[ipart][idim] = 0.0;

  dem->part.tau[ipart] = 0.0;

  dem->cache.stale = true;

  return ipart;
}

/// The call `bmm_dem_copypart(dem, ipart, jpart)`
/// reassigns the particle `jpart` to `ipart`.
__attribute__ ((__nonnull__))
static void bmm_dem_copypart(struct bmm_dem *const dem,
    size_t const ipart, size_t const jpart) {
  dem->part.role[ipart] = dem->part.role[jpart];
  dem->part.l[ipart] = dem->part.l[jpart];
  dem->part.r[ipart] = dem->part.r[jpart];
  dem->part.m[ipart] = dem->part.m[jpart];
  dem->part.jred[ipart] = dem->part.jred[jpart];

  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    dem->part.x[ipart][idim] = dem->part.x[jpart][idim];

  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    dem->part.v[ipart][idim] = dem->part.v[jpart][idim];

  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    dem->part.a[ipart][idim] = dem->part.a[jpart][idim];

  dem->part.phi[ipart] = dem->part.phi[jpart];
  dem->part.omega[ipart] = dem->part.omega[jpart];
  dem->part.alpha[ipart] = dem->part.alpha[jpart];

  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    dem->part.f[ipart][idim] = dem->part.f[jpart][idim];

  dem->part.tau[ipart] = dem->part.tau[jpart];

  for (size_t ict = 0; ict < BMM_NCT; ++ict) {
    dem->pair[ict].cont.src[ipart].n = dem->pair[ict].cont.src[jpart].n;

    for (size_t icont = 0; icont < dem->pair[ict].cont.src[jpart].n; ++icont)
      dem->pair[ict].cont.src[ipart].itgt[icont] = dem->pair[ict].cont.src[jpart].itgt[icont];

    for (size_t icont = 0; icont < dem->pair[ict].cont.src[jpart].n; ++icont)
      dem->pair[ict].cont.src[ipart].drest[icont] = dem->pair[ict].cont.src[jpart].drest[icont];

    for (size_t icont = 0; icont < dem->pair[ict].cont.src[jpart].n; ++icont)
      for (size_t iend = 0; iend < BMM_NEND; ++iend)
        dem->pair[ict].cont.src[ipart].chirest[icont][iend] =
          dem->pair[ict].cont.src[jpart].chirest[icont][iend];

    for (size_t icont = 0; icont < dem->pair[ict].cont.src[jpart].n; ++icont)
      dem->pair[ict].cont.src[ipart].strength[icont] = dem->pair[ict].cont.src[jpart].strength[icont];

    for (size_t icont = 0; icont < dem->pair[ict].cont.src[jpart].n; ++icont)
      dem->pair[ict].cont.src[ipart].rlim[icont] = dem->pair[ict].cont.src[jpart].rlim[icont];

    for (size_t icont = 0; icont < dem->pair[ict].cont.src[jpart].n; ++icont)
      dem->pair[ict].cont.src[ipart].philim[icont] = dem->pair[ict].cont.src[jpart].philim[icont];
  }
}

void bmm_dem_rempart(struct bmm_dem *const dem,
    size_t const ipart) {
  --dem->part.n;

  size_t const jpart = dem->part.n;

  // TODO Go over particles that are in contact with this (N583) and
  // update their associations (otherwise there will be dangling contacts).

  if (jpart != ipart)
    bmm_dem_copypart(dem, ipart, jpart);
}

void bmm_dem_force_creeping(struct bmm_dem *const dem,
    size_t const ipart) {
  double const v = bmm_geom2d_norm(dem->part.v[ipart]);

  if (v == 0.0)
    return;

  double vunit[BMM_NDIM];
  bmm_geom2d_scale(vunit, dem->part.v[ipart], 1.0 / v);

  double const f = -3.0 * M_2PI * dem->amb.params.creeping.mu *
    dem->part.r[ipart] * v;

  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    dem->part.f[ipart][idim] += f * vunit[idim];

  double const tau = -4.0 * M_2PI * dem->amb.params.creeping.mu *
    $(bmm_power, double)(dem->part.r[ipart], 3);

  dem->part.tau[ipart] += tau * dem->part.omega[ipart];
}

// TODO Flatten these to allow loop-invariant code motion.
void bmm_dem_force_ambient(struct bmm_dem *const dem, size_t const ipart) {
  switch (dem->amb.tag) {
    case BMM_DEM_AMB_FAXEN:
      bmm_dem_force_creeping(dem, ipart);

      break;
  }
}

// TODO Don't ever write code like I did here.

void bmm_dem_force_weak(struct bmm_dem *const dem,
    size_t const ipart, size_t const jpart, size_t const icont) {
  double xdiffji[BMM_NDIM];
  bmm_geom2d_cpdiff(xdiffji, dem->part.x[ipart], dem->part.x[jpart],
      dem->opts.box.x, dem->opts.box.per);

  double const d2 = bmm_geom2d_norm2(xdiffji);
  if (d2 == 0.0)
    return;

  double const ri = dem->part.r[ipart];
  double const rj = dem->part.r[jpart];
  double const r = ri + rj;
  double const r2 = $(bmm_power, double)(r, 2);
  if (d2 > r2)
    return;

  double const d = sqrt(d2);

  double xnormji[BMM_NDIM];
  bmm_geom2d_scale(xnormji, xdiffji, 1.0 / d);

  double xtangji[BMM_NDIM];
  bmm_geom2d_rperp(xtangji, xnormji);

  double vdiffij[BMM_NDIM];
  bmm_geom2d_diff(vdiffij, dem->part.v[jpart], dem->part.v[ipart]);

  // Normal forces first.

  double fnorm = 0.0;

  {
    double const xi = r - d;
    double const vnormij = bmm_geom2d_dot(vdiffij, xnormji);
    double const reff = $(bmm_resum2, double)(ri, rj);
    double const dt = dem->opts.script.dt[dem->script.i];

    switch (dem->pair[BMM_DEM_CT_WEAK].norm.tag) {
      case BMM_DEM_NORM_KV:
        fnorm = dem->pair[BMM_DEM_CT_WEAK].norm.params.dashpot.k * xi +
          dem->pair[BMM_DEM_CT_WEAK].norm.params.dashpot.gamma * vnormij;

        break;
      case BMM_DEM_NORM_BSHP:
        fnorm = (2.0 / 3.0) * (dem->opts.part.y /
              (1.0 - $(bmm_power, double)(dem->opts.part.nu, 2))) *
            (xi + dem->pair[BMM_DEM_CT_WEAK].norm.params.viscoel.a * vnormij) * sqrt(reff * xi);

        break;
    }

    fnorm = $(bmm_max, double)(0.0, fnorm);

    // TODO Check this energy.
    dem->est.ewcont += (1.0 / 2.0) * fnorm * vnormij * dt;
  }

  double fnormji[BMM_NDIM];
  bmm_geom2d_scale(fnormji, xnormji, fnorm);

  bmm_geom2d_addto(dem->part.f[ipart], fnormji);
  bmm_geom2d_diffto(dem->part.f[jpart], fnormji);

  // Tangential forces second.

  double ftang = 0.0;

  {
    // TODO By these signs, I think my n--t coordinate system is fucked.
    double const vtangij = bmm_geom2d_dot(vdiffij, xtangji) +
      ri * dem->part.omega[ipart] + rj * dem->part.omega[jpart];
    double const dt = dem->opts.script.dt[dem->script.i];

    switch (dem->pair[BMM_DEM_CT_WEAK].tang.tag) {
      case BMM_DEM_TANG_HW:
        ftang = copysign($(bmm_min, double)(
              dem->pair[BMM_DEM_CT_WEAK].tang.params.hw.gamma * $(bmm_abs, double)(vtangij),
              dem->pair[BMM_DEM_CT_WEAK].tang.params.hw.mu * $(bmm_abs, double)(fnorm)), vtangij);

        break;
      // TODO Not here!
      case BMM_DEM_TANG_CS:
        // TODO Always the case?
        {
          double const epsilon = vtangij * dt;
          if (ipart < jpart)
            dem->pair[BMM_DEM_CT_WEAK].cont.src[ipart].epsilon[icont] += epsilon;
          else
            dem->pair[BMM_DEM_CT_WEAK].cont.src[jpart].epsilon[icont] += epsilon;

          ftang = copysign($(bmm_min, double)(
                dem->pair[BMM_DEM_CT_WEAK].tang.params.cs.kappa * $(bmm_abs, double)(epsilon),
                dem->pair[BMM_DEM_CT_WEAK].tang.params.cs.mu * $(bmm_abs, double)(fnorm)), vtangij);
        }

        break;
    }

    // TODO Check this energy.
    dem->est.ewcont += (1.0 / 2.0) * ftang * vtangij * dt;
  }

  double ftangji[BMM_NDIM];
  bmm_geom2d_scale(ftangji, xtangji, ftang);

  bmm_geom2d_addto(dem->part.f[ipart], ftangji);
  bmm_geom2d_diffto(dem->part.f[jpart], ftangji);

  // Torques third.

  double taui = 0.0;
  double tauj = 0.0;

  {
    // TODO This warrants a sanity check.
    double const ri2 = $(bmm_power, double)(ri, 2);
    double const rj2 = $(bmm_power, double)(rj, 2);
    double const bij = (1.0 / 2.0) * (1.0 - (rj2 - ri2) / d2);
    double const bji = 1.0 - bij;
    double const sij = bij * d;
    double const sji = bji * d;
    double const sij2 = $(bmm_power, double)(sij, 2);
    double const sji2 = $(bmm_power, double)(sji, 2);
    double const cij = 2.0 * sqrt(ri2 - sij2);
    double const cd2sij = cij / (2.0 * sij);
    double const cd2sji = cij / (2.0 * sji);
    double const aij = (1.0 / 2.0) * (ri + (sij / cd2sij) * asinh(cd2sij));
    double const aji = (1.0 / 2.0) * (rj + (sji / cd2sji) * asinh(cd2sji));
    double const aijt = (1.0 / 2.0) * (ri + sij);
    double const ajit = (1.0 / 2.0) * (ri + sji);

    switch (dem->pair[BMM_DEM_CT_WEAK].torque.tag) {
      case BMM_DEM_TORQUE_HARD:
        taui = ri * ftang;
        tauj = rj * ftang;

        break;
      case BMM_DEM_TORQUE_SOFT:
        taui = sij * ftang;
        tauj = sji * ftang;

        break;
      case BMM_DEM_TORQUE_AVERAGE:
        taui = aij * ftang;
        tauj = aji * ftang;

        break;
      case BMM_DEM_TORQUE_HALFWAY:
        taui = aijt * ftang;
        tauj = ajit * ftang;

        break;
    }
  }

  dem->part.tau[ipart] -= taui;
  dem->part.tau[jpart] -= tauj;
}

// TODO This coordinate system has not been flipped yet.

void bmm_dem_force_strong(struct bmm_dem *const dem,
    size_t const ipart, size_t const jpart, size_t const icont) {
  double xdiffij[BMM_NDIM];
  bmm_geom2d_cpdiff(xdiffij, dem->part.x[jpart], dem->part.x[ipart],
      dem->opts.box.x, dem->opts.box.per);

  double const d = bmm_geom2d_norm(xdiffij);
  double const ri = dem->part.r[ipart];
  double const rj = dem->part.r[jpart];
  double const r = ri + rj;
  double const r2 = $(bmm_power, double)(r, 2);

  double xnormij[BMM_NDIM];
  bmm_geom2d_scale(xnormij, xdiffij, 1.0 / d);

  double xtangij[BMM_NDIM];
  bmm_geom2d_rperp(xtangij, xnormij);

  double vdiffij[BMM_NDIM];
  bmm_geom2d_diff(vdiffij, dem->part.v[jpart], dem->part.v[ipart]);

  // Normal forces first.

  double fnorm = 0.0;

  switch (dem->cont_dep.tag) {
    case BMM_DEM_LINK_SPRING:
    case BMM_DEM_LINK_BEAM:
      {
        double const d0 = dem->pair[BMM_DEM_CT_STRONG].cont.src[ipart].drest[icont];
        double const drest = d * dem->opts.cont.cshcont;
        double const zeta = d0 - d;
        double const vnormji = -bmm_geom2d_dot(vdiffij, xnormij);

        fnorm = -(dem->opts.cont.ktens * zeta + dem->opts.cont.dktens * vnormji);

        double fnormij[BMM_NDIM];
        bmm_geom2d_scale(fnormij, xnormij, fnorm);

        double fnormji[BMM_NDIM];
        bmm_geom2d_scale(fnormji, fnormij, -1.0);

        bmm_geom2d_addto(dem->part.f[ipart], fnormij);
        bmm_geom2d_addto(dem->part.f[jpart], fnormji);

        // TODO Check this energy.
        dem->est.escont += (1.0 / 2.0) * dem->opts.cont.ktens *
            $(bmm_power, double)(zeta / d0, 2);
      }
      break;
  }

  switch (dem->cont_dep.tag) {
    case BMM_DEM_LINK_BEAM:
      {
        // Torques second.

        double taui = 0.0;
        double tauj = 0.0;

        double const gammaij = $(bmm_swrap, double)(bmm_geom2d_dir(xdiffij), M_2PI);
        double const gammaji = $(bmm_swrap, double)(gammaij + M_PI, M_2PI);
        double const phii = $(bmm_swrap, double)(dem->part.phi[ipart], M_2PI);
        double const phij = $(bmm_swrap, double)(dem->part.phi[jpart], M_2PI);
        double const chii = $(bmm_swrap, double)(gammaij - phii, M_2PI);
        double const chij = $(bmm_swrap, double)(gammaji - phij, M_2PI);

        double const dchii = $(bmm_swrap, double)(dem->pair[BMM_DEM_CT_STRONG].cont.src[ipart].chirest[icont][BMM_DEM_END_TAIL] - chii, M_2PI);
        double const dchij = $(bmm_swrap, double)(dem->pair[BMM_DEM_CT_STRONG].cont.src[ipart].chirest[icont][BMM_DEM_END_HEAD] - chij, M_2PI);

        // TODO Derive these on paper too.
        double const vrotij = -bmm_geom2d_dot(vdiffij, xtangij);
        double const omegai = dem->part.omega[ipart] + vrotij / d;
        double const omegaj = dem->part.omega[jpart] + vrotij / d;

        taui = -(dem->opts.cont.kshear * dchii + dem->opts.cont.dkshear * omegai);
        tauj = -(dem->opts.cont.kshear * dchij + dem->opts.cont.dkshear * omegaj);

        dem->part.tau[ipart] += taui;
        dem->part.tau[jpart] += tauj;

        // Tangential forces third.

        double ftang = 0.0;

        {
          ftang = (taui + tauj) / d;
        }

        double ftangij[BMM_NDIM];
        bmm_geom2d_scale(ftangij, xtangij, ftang);

        double ftangji[BMM_NDIM];
        bmm_geom2d_scale(ftangji, ftangij, -1.0);

        bmm_geom2d_addto(dem->part.f[ipart], ftangij);
        bmm_geom2d_addto(dem->part.f[jpart], ftangji);

        // TODO Check this energy.
        dem->est.escont += (1.0 / 2.0) * dem->opts.cont.kshear *
            $(bmm_power, double)(dchij, 2);
      }

      break;
  }
}

void bmm_dem_force_external(struct bmm_dem *const dem, size_t const ipart) {
  switch (dem->ext.tag) {
    case BMM_DEM_EXT_HARM:
      dem->part.f[ipart][1] += dem->ext.params.harm.k *
        (dem->opts.box.x[1] / 2.0 - dem->part.x[ipart][1]);

      break;
    case BMM_DEM_EXT_GRAVY:
      dem->part.f[ipart][1] += dem->part.m[ipart] * dem->ext.params.gravy.g;

      break;
    case BMM_DEM_EXT_DRIVE:
      if (dem->part.role[ipart] == BMM_DEM_ROLE_DRIVEN)
        for (size_t idim = 0; idim < BMM_NDIM; ++idim)
          dem->part.f[ipart][idim] +=
            dem->script.state.crunch.f[idim] /
            (double) dem->script.state.crunch.ndrive;

      break;
  }
}

void bmm_dem_force(struct bmm_dem *const dem) {
  double const dt = dem->opts.script.dt[dem->script.i];

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    for (size_t idim = 0; idim < BMM_NDIM; ++idim)
      dem->part.f[ipart][idim] = 0.0;

    dem->part.tau[ipart] = 0.0;
  }

  // Estimators appeared.
  dem->est.epotext = NAN;
  dem->est.eklin = NAN;
  dem->est.ekrot = NAN;
  dem->est.ewcont = 0.0;
  dem->est.escont = 0.0;
  dem->est.edrivnorm = 0.0;
  dem->est.edrivtang = 0.0;
  dem->est.eyieldis = 0.0;
  dem->est.econtdis = 0.0;
  dem->est.efricdis = 0.0;
  dem->est.fback[0] = 0.0;
  dem->est.fback[1] = 0.0;
  dem->est.mueff = 0.0;
  dem->est.vdriv[0] = NAN;
  dem->est.vdriv[1] = NAN;

  // This goes for the previous frame, so this is not exactly the right spot.
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    if (dem->part.role[ipart] == BMM_DEM_ROLE_DRIVEN)
      for (size_t idim = 0; idim < BMM_NDIM; ++idim)
        dem->est.vdriv[idim] += dem->part.v[ipart][idim];

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    bmm_dem_force_ambient(dem, ipart);

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    for (size_t icont = 0; icont < dem->pair[BMM_DEM_CT_STRONG].cont.src[ipart].n; ++icont) {
      size_t const jpart = dem->pair[BMM_DEM_CT_STRONG].cont.src[ipart].itgt[icont];

      bmm_dem_force_strong(dem, ipart, jpart, icont);
    }

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    for (size_t icont = 0; icont < dem->pair[BMM_DEM_CT_WEAK].cont.src[ipart].n; ++icont) {
      size_t const jpart = dem->pair[BMM_DEM_CT_WEAK].cont.src[ipart].itgt[icont];

      bmm_dem_force_weak(dem, ipart, jpart, icont);
    }

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    if (dem->part.role[ipart] == BMM_DEM_ROLE_DRIVEN)
      for (size_t idim = 0; idim < BMM_NDIM; ++idim)
        dem->est.fback[idim] += dem->part.f[ipart][idim];

  if (dem->est.fback[1] != 0.0)
    dem->est.mueff = $(bmm_abs, double)(dem->est.fback[0] / dem->est.fback[1]);

  double v[BMM_NDIM];
  bmm_dem_est_vdc(dem->est.vdriv, dem);

  dem->script.state.crunch.f[0] +=
    copysign(dem->opts.script.params[dem->script.i].crunch.fadjust[0] * dt,
        dem->opts.script.params[dem->script.i].crunch.v - dem->est.vdriv[0]);

  // TODO Not like this.
  dem->script.state.crunch.f[1] = dem->opts.script.params[dem->script.i].crunch.p;

  // TODO This still sucks.
  switch (dem->ext.tag) {
    case BMM_DEM_EXT_DRIVE:
      dem->script.state.crunch.ndrive = 0;
      for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
        if (dem->part.role[ipart] == BMM_DEM_ROLE_DRIVEN)
          ++dem->script.state.crunch.ndrive;

      break;
  }
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    bmm_dem_force_external(dem, ipart);
}

void bmm_dem_accel(struct bmm_dem *const dem) {
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    double const m = dem->part.m[ipart];
    double const j = dem->cache.j[ipart];

    for (size_t idim = 0; idim < BMM_NDIM; ++idim)
      dem->part.a[ipart][idim] = dem->part.f[ipart][idim] / m;

    dem->part.alpha[ipart] = dem->part.tau[ipart] / j;
  }
}

void bmm_dem_integ_euler(struct bmm_dem *const dem) {
  double const dt = dem->opts.script.dt[dem->script.i];

  if (dt == 0.0)
    return;

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    if (dem->part.role[ipart] != BMM_DEM_ROLE_FIXED) {
      for (size_t idim = 0; idim < BMM_NDIM; ++idim) {
        dem->part.x[ipart][idim] +=
          (1.0 / (double) BMM_FACT(1)) * dem->part.v[ipart][idim] * dt;
        dem->part.v[ipart][idim] +=
          (1.0 / (double) BMM_FACT(1)) * dem->part.a[ipart][idim] * dt;

        if (dem->opts.box.per[idim])
          dem->part.x[ipart][idim] = $(bmm_uwrap, double)(dem->part.x[ipart][idim],
              dem->opts.box.x[idim]);
      }

      dem->part.phi[ipart] +=
        (1.0 / (double) BMM_FACT(1)) * dem->part.omega[ipart] * dt;
      dem->part.omega[ipart] +=
        (1.0 / (double) BMM_FACT(1)) * dem->part.alpha[ipart] * dt;
    }
}

void bmm_dem_integ_taylor(struct bmm_dem *const dem) {
  double const dt = dem->opts.script.dt[dem->script.i];

  if (dt == 0.0)
    return;

  double const dt2 = $(bmm_power, double)(dt, 2);

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    if (dem->part.role[ipart] != BMM_DEM_ROLE_FIXED) {
      for (size_t idim = 0; idim < BMM_NDIM; ++idim) {
        dem->part.x[ipart][idim] +=
          (1.0 / (double) BMM_FACT(1)) * dem->part.v[ipart][idim] * dt +
          (1.0 / (double) BMM_FACT(2)) * dem->part.a[ipart][idim] * dt2;
        dem->part.v[ipart][idim] +=
          (1.0 / (double) BMM_FACT(1)) * dem->part.a[ipart][idim] * dt;

        if (dem->opts.box.per[idim])
          dem->part.x[ipart][idim] = $(bmm_uwrap, double)(dem->part.x[ipart][idim],
              dem->opts.box.x[idim]);
      }

      dem->part.phi[ipart] +=
        (1.0 / (double) BMM_FACT(1)) * dem->part.omega[ipart] * dt +
        (1.0 / (double) BMM_FACT(2)) * dem->part.alpha[ipart] * dt2;
      dem->part.omega[ipart] +=
        (1.0 / (double) BMM_FACT(1)) * dem->part.alpha[ipart] * dt;
    }
}

void bmm_dem_integ_vel(struct bmm_dem *const dem) {
  double const dt = dem->opts.script.dt[dem->script.i];

  if (dt == 0.0)
    return;

  double const dt2 = $(bmm_power, double)(dt, 2);

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    if (dem->part.role[ipart] != BMM_DEM_ROLE_FIXED) {
      for (size_t idim = 0; idim < BMM_NDIM; ++idim) {
        dem->part.x[ipart][idim] +=
          (1.0 / (double) BMM_FACT(1)) * dem->part.v[ipart][idim] * dt +
          (1.0 / (double) BMM_FACT(2)) * dem->part.a[ipart][idim] * dt2;

        if (dem->opts.box.per[idim])
          dem->part.x[ipart][idim] = $(bmm_uwrap, double)(dem->part.x[ipart][idim],
              dem->opts.box.x[idim]);

        dem->integ.params.velvet.a[ipart][idim] = dem->part.a[ipart][idim];
      }

      dem->part.phi[ipart] +=
        (1.0 / (double) BMM_FACT(1)) * dem->part.omega[ipart] * dt +
        (1.0 / (double) BMM_FACT(2)) * dem->part.alpha[ipart] * dt2;

      dem->integ.params.velvet.alpha[ipart] = dem->part.alpha[ipart];
    }
}

void bmm_dem_integ_vet(struct bmm_dem *const dem) {
  double const dt = dem->opts.script.dt[dem->script.i];

  if (dt == 0.0)
    return;

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    if (dem->part.role[ipart] != BMM_DEM_ROLE_FIXED) {
      for (size_t idim = 0; idim < BMM_NDIM; ++idim)
        dem->part.v[ipart][idim] += (1.0 / (double) BMM_FACT(2)) *
          (dem->integ.params.velvet.a[ipart][idim] + dem->part.a[ipart][idim]) * dt;

      dem->part.omega[ipart] += (1.0 / (double) BMM_FACT(2)) *
        (dem->integ.params.velvet.alpha[ipart] + dem->part.alpha[ipart]) * dt;
    }
}

// TODO This erases the winding number information needed by beams.
void bmm_dem_stab(struct bmm_dem *const dem) {
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    dem->part.phi[ipart] = $(bmm_uwrap, double)(dem->part.phi[ipart], M_2PI);
}

void bmm_dem_predict(struct bmm_dem *const dem) {
  switch (dem->integ.tag) {
    case BMM_DEM_INTEG_VELVET:
      bmm_dem_integ_vel(dem);

      break;
  }
}

void bmm_dem_correct(struct bmm_dem *const dem) {
  switch (dem->integ.tag) {
    case BMM_DEM_INTEG_EULER:
      bmm_dem_integ_euler(dem);

      break;
    case BMM_DEM_INTEG_TAYLOR:
      bmm_dem_integ_taylor(dem);

      break;
    case BMM_DEM_INTEG_VELVET:
      bmm_dem_integ_vet(dem);

      break;
  }
}

bool bmm_dem_fault_pair(struct bmm_dem *const dem,
    size_t const ipart, size_t const jpart, size_t const icont) {
  bool const iind = dem->opts.script.params[dem->script.i].fault.
    ind(dem->part.x[ipart], &dem->opts);

  bool const jind = dem->opts.script.params[dem->script.i].fault.
    ind(dem->part.x[jpart], &dem->opts);

  if (iind != jind) {
    bmm_dem_remcont(dem, BMM_DEM_CT_STRONG, ipart, jpart, icont);

    return true;
  }

  return false;
}

void bmm_dem_fault(struct bmm_dem *const dem) {
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    for (size_t icont = 0; icont < dem->pair[BMM_DEM_CT_STRONG].cont.src[ipart].n; ++icont) {
      size_t const jpart = dem->pair[BMM_DEM_CT_STRONG].cont.src[ipart].itgt[icont];

      if (bmm_dem_fault_pair(dem, ipart, jpart, icont))
        --icont;
    }
}

bool bmm_dem_link_pair(struct bmm_dem *const dem,
    size_t const ipart, size_t const jpart) {
  // TODO Settle these order problems.
  if (ipart >= jpart) {
    return bmm_dem_link_pair(dem, jpart, ipart);
  }

  double xdiffij[BMM_NDIM];
  bmm_geom2d_cpdiff(xdiffij, dem->part.x[jpart], dem->part.x[ipart],
      dem->opts.box.x, dem->opts.box.per);

  double const d2 = bmm_geom2d_norm2(xdiffij);

  double const r = dem->part.r[ipart] + dem->part.r[jpart];
  double const r2 = $(bmm_power, double)(r, 2);

  if (d2 > r2 * dem->opts.cont.ccrcont)
    return true;

  size_t const icont = bmm_dem_search_cont(dem, BMM_DEM_CT_WEAK, ipart, jpart);
  if (icont != SIZE_MAX)
    bmm_dem_remcont(dem, BMM_DEM_CT_WEAK, ipart, jpart, icont);

  return bmm_dem_addcont(dem, BMM_DEM_CT_STRONG, ipart, jpart) != SIZE_MAX;
}

// TODO Check triangulation quality and compare with Delaunay.
bool bmm_dem_link(struct bmm_dem *const dem) {
  switch (dem->cache.tag) {
    case BMM_DEM_CACHE_NONE:
      for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
        for (size_t jpart = ipart + 1; jpart < dem->part.n; ++jpart)
          if (!bmm_dem_link_pair(dem, ipart, jpart))
            return false;

      break;
    case BMM_DEM_CACHE_NEIGH:
      for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
        for (size_t ineigh = 0; ineigh < dem->cache.neigh[ipart].n; ++ineigh)
          if (!bmm_dem_link_pair(dem, ipart, dem->cache.neigh[ipart].i[ineigh]))
            return false;

      break;
  }

  return true;
}

bool bmm_dem_uncont(struct bmm_dem *const dem) {
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    dem->pair[BMM_DEM_CT_STRONG].cont.src[ipart].n = 0;

  return true;
}

// TODO Interval arithmetic.
// TODO Not quite! Edges may overlap if periodicity is turned on.
__attribute__ ((__nonnull__, __pure__))
static bool bmm_dem_inside(struct bmm_dem const *const dem,
    size_t const ipart) {
  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    if (dem->part.x[ipart][idim] < 0.0 ||
        dem->part.x[ipart][idim] >= dem->opts.box.x[idim])
      return false;

  return true;
}

// TODO Not quite! Edges may be excavated if periodicity is turned off.
__attribute__ ((__nonnull__, __pure__))
static bool bmm_dem_insider(struct bmm_dem const *const dem,
    size_t const ipart) {
  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    if (dem->part.x[ipart][idim] - dem->part.r[ipart] < 0.0 ||
        dem->part.x[ipart][idim] + dem->part.r[ipart] >= dem->opts.box.x[idim])
      return false;

  return true;
}

__attribute__ ((__nonnull__))
static void bmm_dem_script_clip(struct bmm_dem *const dem) {
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    if (!bmm_dem_inside(dem, ipart)) {
      bmm_dem_rempart(dem, ipart);
      --ipart;
    }

  // TODO See (N583).
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    for (size_t ict = 0; ict < BMM_NCT; ++ict)
      dem->pair[ict].cont.src[ipart].n = 0;
}

static bool extra_crap(struct bmm_dem const *const dem) {
  double const epotext = dem->est.epotext;
  double const eklin = dem->est.eklin;
  double const ekrot = dem->est.ekrot;
  double const ewcont = dem->est.ewcont;
  double const escont = dem->est.escont;
  double const edrivnorm = dem->est.edrivnorm;
  double const edrivtang = dem->est.edrivtang;
  double const eyieldis = dem->est.eyieldis;
  double const econtdis = dem->est.econtdis;
  double const efricdis = dem->est.efricdis;
  double const pos = epotext + eklin + ekrot + ewcont + escont + edrivnorm + edrivtang;
  double const neg = eyieldis + econtdis + efricdis;
  if (fprintf(stderr,
        "epotext = %g, "
        "eklin = %g, "
        "ekrot = %g, "
        "ewcont = %g, "
        "escont = %g, "
        "edrivnorm = %g, "
        "edrivtang = %g, "
        "eyieldis = %g, "
        "econtdis = %g, "
        "efricdis = %g, "
        "sum(+) = %g, "
        "sum(-) = %g, "
        "force = %g, "
        "sum(?) = %g\n",
        epotext, eklin, ekrot, ewcont, escont, edrivnorm, edrivtang,
        eyieldis, econtdis, efricdis,
        pos, neg,
        dem->script.state.crunch.f[0],
        pos - neg) < 0) {
    BMM_TLE_STDS();

    return false;
  }

  return true;
}

static bool export_x(struct bmm_dem const *const dem) {
  FILE *const stream = fopen("fragment-x.data", "w");
  if (stream == NULL) {
    BMM_TLE_STDS();

    return false;
  }

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    if (fprintf(stream, "(%g, %g)\n",
          dem->part.x[ipart][0], dem->part.x[ipart][1]) < 0) {
      BMM_TLE_STDS();

      return false;
    }

  if (fclose(stream) != 0) {
    BMM_TLE_STDS();

    return false;
  }

  return true;
}

static bool export_r(struct bmm_dem const *const dem) {
  FILE *const stream = fopen("fragment-r.data", "w");
  if (stream == NULL) {
    BMM_TLE_STDS();

    return false;
  }

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    if (fprintf(stream, "%g\n", dem->part.r[ipart]) < 0) {
      BMM_TLE_STDS();

      return false;
    }

  if (fclose(stream) != 0) {
    BMM_TLE_STDS();

    return false;
  }

  return true;
}

static bool export_c(struct bmm_dem const *const dem) {
  FILE *const stream = fopen("fragment-c.data", "w");
  if (stream == NULL) {
    BMM_TLE_STDS();

    return false;
  }

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    if (fprintf(stream, "%d\n", dem->part.role[ipart]) < 0) {
      BMM_TLE_STDS();

      return false;
    }

  if (fclose(stream) != 0) {
    BMM_TLE_STDS();

    return false;
  }

  return true;
}

struct agraph {
  struct {
    size_t n;
    size_t itgt[BMM_MCONTACT * 2];
    bool visited[BMM_MCONTACT * 2];
  } src[BMM_MPART];
};

struct acls {
  struct agraph *agraph;
  struct bmm_dem const *dem;
  size_t ipart;
};

struct aface {
  struct {
    size_t n;
    size_t ivert[BMM_MPART * 2];
    double circ;
    double area;
  } poly[BMM_MPART * 4];
  size_t lastseen[BMM_MPART];
};

__attribute__ ((__nonnull__, __pure__))
static int acompar(size_t const i, size_t const j, void *const cls) {
  struct acls const *const acls = cls;
  struct bmm_dem const *const dem = acls->dem;
  struct agraph const *const agraph = acls->agraph;
  size_t const ipart = acls->ipart;

  double gammai;
  {
    size_t const jpart = agraph->src[ipart].itgt[i];

    double xdiffij[BMM_NDIM];
    bmm_geom2d_diff(xdiffij, dem->part.x[jpart], dem->part.x[ipart]);

    gammai = $(bmm_swrap, double)(bmm_geom2d_dir(xdiffij), M_2PI);
  }

  double gammaj;
  {
    size_t const jpart = agraph->src[ipart].itgt[j];

    double xdiffij[BMM_NDIM];
    bmm_geom2d_diff(xdiffij, dem->part.x[jpart], dem->part.x[ipart]);

    gammaj = $(bmm_swrap, double)(bmm_geom2d_dir(xdiffij), M_2PI);
  }

  return $(bmm_cmp, double)(gammaj, gammai);
}

__attribute__ ((__nonnull__))
static void aswap(size_t const i, size_t const j, void *const cls) {
  struct acls *const acls = cls;
  struct bmm_dem const *const dem = acls->dem;
  struct agraph *const agraph = acls->agraph;
  size_t const ipart = acls->ipart;

  size_t const tmp = agraph->src[ipart].itgt[i];
  agraph->src[ipart].itgt[i] = agraph->src[ipart].itgt[j];
  agraph->src[ipart].itgt[j] = tmp;
}

static bool export_f(struct bmm_dem const *const dem) {
  FILE *const stream = fopen("fragment-f.data", "w");
  if (stream == NULL) {
    BMM_TLE_STDS();

    return false;
  }

  // Bidirectionalization of the strong contact graph.
  // Following great programming conventions, I copied this from below.
  struct agraph *const agraph = malloc(sizeof *agraph);
  if (agraph == NULL)
    BMM_TLE_STDS();
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    agraph->src[ipart].n = 0;
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    for (size_t icont = 0; icont < dem->pair[BMM_DEM_CT_STRONG].cont.src[ipart].n; ++icont) {
      size_t const jpart = dem->pair[BMM_DEM_CT_STRONG].cont.src[ipart].itgt[icont];

      double d[BMM_NDIM];
      for (size_t idim = 0; idim < BMM_NDIM; ++idim)
        d[idim] = $(bmm_abs, double)
          (dem->part.x[jpart][idim] - dem->part.x[ipart][idim]);

      // No winding.
      if (d[0] < dem->opts.box.x[0] / 2.0 &&
          d[1] < dem->opts.box.x[1] / 2.0) {
        size_t const jcont = agraph->src[jpart].n;
        agraph->src[jpart].itgt[jcont] = ipart;
        ++agraph->src[jpart].n;

        // This makes them bidirectional.
        size_t const kcont = agraph->src[ipart].n;
        agraph->src[ipart].itgt[kcont] = jpart;
        ++agraph->src[ipart].n;
      }
    }

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    if (agraph->src[ipart].n == 0)
      if (fprintf(stream, "%zu\n", ipart) < 0) {
        BMM_TLE_STDS();

        return false;
      }

  free(agraph);

  if (fclose(stream) != 0) {
    BMM_TLE_STDS();

    return false;
  }

  return true;
}

static bool export_p(struct bmm_dem const *const dem) {
  FILE *const stream = fopen("fragment-p.data", "w");
  if (stream == NULL) {
    BMM_TLE_STDS();

    return false;
  }

  // Bidirectionalization of the strong contact graph.
  struct agraph *const agraph = malloc(sizeof *agraph);
  if (agraph == NULL)
    BMM_TLE_STDS();
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    agraph->src[ipart].n = 0;
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    for (size_t icont = 0; icont < dem->pair[BMM_DEM_CT_STRONG].cont.src[ipart].n; ++icont) {
      size_t const jpart = dem->pair[BMM_DEM_CT_STRONG].cont.src[ipart].itgt[icont];

      double d[BMM_NDIM];
      for (size_t idim = 0; idim < BMM_NDIM; ++idim)
        d[idim] = $(bmm_abs, double)
          (dem->part.x[jpart][idim] - dem->part.x[ipart][idim]);

      // No winding.
      if (d[0] < dem->opts.box.x[0] / 2.0 &&
          d[1] < dem->opts.box.x[1] / 2.0) {
        size_t const jcont = agraph->src[jpart].n;
        agraph->src[jpart].itgt[jcont] = ipart;
        ++agraph->src[jpart].n;

        // This makes them bidirectional.
        size_t const kcont = agraph->src[ipart].n;
        agraph->src[ipart].itgt[kcont] = jpart;
        ++agraph->src[ipart].n;
      }
    }

  // Angle ordering.
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    struct acls acls = {
      .agraph = agraph,
      .dem = dem,
      .ipart = ipart
    };
    bmm_hsort_cls(agraph->src[ipart].n, acompar, aswap, &acls);
    /*
    for (size_t iarr = 0; iarr < agraph->src[ipart].n / 2; ++iarr)
      agraph->src[ipart].itgt[iarr] = agraph->src[ipart].itgt[agraph->src[ipart].n - 1 - iarr];
    */
  }

  // Preparations.
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    for (size_t iarr = 0; iarr < agraph->src[ipart].n; ++iarr)
      agraph->src[ipart].visited[iarr] = false;

  // Elimination walk.
  struct aface *const faces = malloc(sizeof *faces);
  if (faces == NULL)
    BMM_TLE_STDS();
  size_t nface = 0;
  size_t ipart = 0;
  while (ipart < dem->part.n) {
    if (nface >= nmembof(faces->poly)) {
      fprintf(stderr, "Out of face.\n");
      break;
    }

    faces->poly[nface].n = 0;

    faces->poly[nface].ivert[faces->poly[nface].n] = ipart;
    ++faces->poly[nface].n;

    size_t narr = 0;
    for (size_t iarr = 0; iarr < agraph->src[ipart].n; ++iarr)
      if (!agraph->src[ipart].visited[iarr])
        ++narr;
    if (narr == 0) {
      ++nface;

      ++ipart;

      continue;
    }

    size_t prev = SIZE_MAX;
    size_t kpart = ipart;

    size_t initarr = 0;
    for (size_t iarr = 0; iarr < agraph->src[kpart].n; ++iarr) {
      if (!agraph->src[kpart].visited[iarr]) {
        initarr = iarr;

        break;
      }
    }

more: ;

    size_t inext;
    if (prev == SIZE_MAX)
      inext = initarr;
    else {
      size_t iback = 0;

      for (size_t iarr = 0; iarr < agraph->src[kpart].n; ++iarr)
        if (agraph->src[kpart].itgt[iarr] == prev) {
          iback = iarr;

          break;
        }

      inext = 0;
      for (size_t iarr = 0; iarr < agraph->src[kpart].n; ++iarr) {
        size_t const off = (iback + 1 + iarr) % agraph->src[kpart].n;
        if (!agraph->src[kpart].visited[off]) {
          inext = off;

          break;
        }
      }
    }

    // fprintf(stderr, "Visited %zu -> %zu.\n", kpart, agraph->src[kpart].itgt[inext]);
    size_t const jpart = agraph->src[kpart].itgt[inext];

    if (prev != SIZE_MAX)
      agraph->src[kpart].visited[inext] = true;

    if (jpart != ipart) {
      prev = kpart;
      kpart = jpart;

      faces->poly[nface].ivert[faces->poly[nface].n] = kpart;
      ++faces->poly[nface].n;

      goto more;
    } else {
      size_t inextnext;
      if (kpart == SIZE_MAX)
        inextnext = initarr;
      else {
        size_t ibackback = 0;

        for (size_t iarr = 0; iarr < agraph->src[jpart].n; ++iarr)
          if (agraph->src[jpart].itgt[iarr] == kpart) {
            ibackback = iarr;

            break;
          }

        inextnext = 0;
        for (size_t iarr = 0; iarr < agraph->src[jpart].n; ++iarr) {
          size_t const off = (ibackback + 1 + iarr) % agraph->src[jpart].n;
          if (!agraph->src[jpart].visited[off]) {
            inextnext = off;

            break;
          }
        }
      }

      if (inextnext != initarr) {
        prev = kpart;
        kpart = jpart;

        faces->poly[nface].ivert[faces->poly[nface].n] = kpart;
        ++faces->poly[nface].n;

        goto more;
      } else {
        agraph->src[ipart].visited[initarr] = true;

        ++nface;

        continue;
      }
    }
  }

  // Split self-intersections.
restart: ;
  for (size_t iface = 0; iface < nface; ++iface) {
    for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
      faces->lastseen[ipart] = SIZE_MAX;

    /*
    fprintf(stderr, "Before split face %zu:", iface);
    for (size_t ivert = 0; ivert < faces->poly[iface].n; ++ivert)
      fprintf(stderr, " %zu", faces->poly[iface].ivert[ivert]);
    fprintf(stderr, ".\n");
    */

    for (size_t ivert = 0; ivert < faces->poly[iface].n; ++ivert) {
      size_t const ipart = faces->poly[iface].ivert[ivert];
      size_t const jvert = faces->lastseen[ipart];

      if (jvert == ivert) {
        fprintf(stderr, "Choked real bad.\n");
        goto unstart;
      }

      if (jvert != SIZE_MAX) {
        /*
        fprintf(stderr, "Needs to split at %zu (%zu) to %zu (%zu).\n",
            jvert, faces->poly[iface].ivert[jvert],
            ivert, faces->poly[iface].ivert[ivert]);
        */

        if (nface >= nmembof(faces->poly)) {
          fprintf(stderr, "Out of face again.\n");
          goto unstart;
        }

        // Copy midsection to new head.
        faces->poly[nface].n = 0;
        for (size_t iarr = jvert; iarr < ivert; ++iarr) {
          faces->poly[nface].ivert[faces->poly[nface].n] = faces->poly[iface].ivert[iarr];
          ++faces->poly[nface].n;
        }
        ++nface;

        // Move tail over midsection.
        for (size_t iarr = ivert; iarr < faces->poly[iface].n; ++iarr)
          faces->poly[iface].ivert[iarr - (ivert - jvert)] = faces->poly[iface].ivert[iarr];
        faces->poly[iface].n = faces->poly[iface].n - (ivert - jvert);

        /*
        fprintf(stderr, "After split face %zu: ", iface);
        for (size_t ivert = 0; ivert < faces->poly[iface].n; ++ivert)
          fprintf(stderr, "%zu ", faces->poly[iface].ivert[ivert]);
        fprintf(stderr, "end.\n");

        fprintf(stderr, "After split face %zu: ", nface - 1);
        for (size_t ivert = 0; ivert < faces->poly[nface - 1].n; ++ivert)
          fprintf(stderr, "%zu ", faces->poly[nface - 1].ivert[ivert]);
        fprintf(stderr, "end.\n");
        */

        goto restart;
      } else
        faces->lastseen[ipart] = ivert;
    }
  }
unstart: ;

  // Face metrics.
  for (size_t iface = 0; iface < nface; ++iface) {
    faces->poly[iface].circ = 0.0;
    faces->poly[iface].area = 0.0;

    for (size_t ivert = 1; ivert < faces->poly[iface].n; ++ivert) {
      size_t const ipart = faces->poly[iface].ivert[ivert - 1];
      size_t const jpart = faces->poly[iface].ivert[ivert];

      double xdiffij[BMM_NDIM];
      bmm_geom2d_diff(xdiffij, dem->part.x[jpart], dem->part.x[ipart]);

      faces->poly[iface].circ += bmm_geom2d_norm(xdiffij);
    }

    for (size_t ilace = 0; ilace < faces->poly[iface].n; ++ilace) {
      size_t const ivert = ilace % faces->poly[iface].n;
      size_t const jvert = (ilace + 1) % faces->poly[iface].n;

      size_t const ipart = faces->poly[iface].ivert[ivert];
      size_t const jpart = faces->poly[iface].ivert[jvert];
      // fprintf(stderr, "Trilaterated %zu -> %zu.\n", ipart, jpart);

      faces->poly[iface].area +=
        dem->part.x[ipart][0] * dem->part.x[jpart][1] -
        dem->part.x[jpart][0] * dem->part.x[ipart][1];
    }

    faces->poly[iface].area /= 2.0;
    // fprintf(stderr, "Face %zu with area %g.\n", iface, faces->poly[iface].area);
  }

  // Filter.
  for (size_t iface = 0; iface < nface; ++iface)
    // This would be a useful feature,
    // but with the current algorithm it loses edges
    // that are surrounded by a hole and the outside.
    /*
    if (faces->poly[iface].area < -1.0e-12 ||
        faces->poly[iface].circ > 8.0 * bmm_ival_midpoint(dem->opts.part.rnew)) {
    */
    if (faces->poly[iface].area < -1.0e-12) {
      // fprintf(stderr, "Removed face %zu.\n", iface);
      --nface;
      faces->poly[iface].n = faces->poly[nface].n;
      faces->poly[iface].circ = faces->poly[nface].circ;
      faces->poly[iface].area = faces->poly[nface].area;
      for (size_t ivert = 0; ivert < faces->poly[nface].n; ++ivert)
        faces->poly[iface].ivert[ivert] = faces->poly[nface].ivert[ivert];

      --iface;
    }

  // Rebuild.
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    agraph->src[ipart].n = 0;
  for (size_t iface = 0; iface < nface; ++iface)
    for (size_t ivert = 0; ivert < faces->poly[iface].n; ++ivert) {
      size_t const jvert = (ivert + 1) % faces->poly[iface].n;
      size_t const ipart = faces->poly[iface].ivert[ivert];
      size_t const jpart = faces->poly[iface].ivert[jvert];

      if (jpart != ipart) {
        size_t const icont = agraph->src[ipart].n;
        agraph->src[ipart].itgt[icont] = jpart;
        ++agraph->src[ipart].n;
      }
    }

  // Preparations.
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    for (size_t iarr = 0; iarr < agraph->src[ipart].n; ++iarr)
      agraph->src[ipart].visited[iarr] = false;

  // Double edge removal.
harder: ;
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    for (size_t iarr = 0; iarr < agraph->src[ipart].n; ++iarr) {
      size_t const jpart = agraph->src[ipart].itgt[iarr];

      for (size_t jarr = 0; jarr < agraph->src[jpart].n; ++jarr)
        if (agraph->src[jpart].itgt[jarr] == ipart) {
          if (!agraph->src[ipart].visited[iarr] &&
              !agraph->src[jpart].visited[jarr]) {
            // fprintf(stderr, "Double removed %zu -> %zu.\n", ipart, jpart);
            agraph->src[ipart].visited[iarr] = true;
            agraph->src[jpart].visited[jarr] = true;

            goto harder;
          }
        }
    }

  // Restore zero-area double edges.
  for (size_t iface = 0; iface < nface; ++iface)
    if (fabs(faces->poly[iface].area) < 1.0e-12) {
      for (size_t ivert = 0; ivert < faces->poly[iface].n; ++ivert) {
        size_t const jvert = (ivert + 1) % faces->poly[iface].n;
        size_t const ipart = faces->poly[iface].ivert[ivert];
        size_t const jpart = faces->poly[iface].ivert[jvert];

        for (size_t iarr = 0; iarr < agraph->src[ipart].n; ++iarr)
          if (agraph->src[ipart].itgt[iarr] == jpart) {
            for (size_t jarr = 0; jarr < agraph->src[jpart].n; ++jarr)
              if (agraph->src[jpart].itgt[jarr] == ipart) {
                // fprintf(stderr, "Double restored %zu -> %zu.\n", ipart, jpart);

                agraph->src[ipart].visited[iarr] = false;
                agraph->src[jpart].visited[jarr] = false;

                break;
              }

            break;
          }
      }
    }

  /*
  fprintf(stderr, "digraph {\n");
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    for (size_t iarr = 0; iarr < agraph->src[ipart].n; ++iarr) {
      size_t const jpart = agraph->src[ipart].itgt[iarr];

      if (!agraph->src[ipart].visited[iarr])
        fprintf(stderr, "%zu -> %zu\n", ipart, jpart);
    }
  fprintf(stderr, "}\n");
  */

  // Connected components.
  nface = 0;
  ipart = 0;
  while (ipart < dem->part.n) {
    if (nface >= nmembof(faces->poly)) {
      fprintf(stderr, "Out of face.\n");
      break;
    }

    faces->poly[nface].n = 0;

    double backdir = NAN;
    size_t kpart = ipart;

morer: ;

    size_t narr = 0;
    for (size_t iarr = 0; iarr < agraph->src[kpart].n; ++iarr)
      if (!agraph->src[kpart].visited[iarr])
        ++narr;
    if (narr == 0) {
      // fprintf(stderr, "Face %zu done! ", nface);

      if (faces->poly[nface].n != 0)
        ++nface;

      ++ipart;

      continue;
    }

    faces->poly[nface].ivert[faces->poly[nface].n] = kpart;
    ++faces->poly[nface].n;

    size_t inext = SIZE_MAX;
    if (isnan(backdir)) {
      for (size_t iarr = 0; iarr < agraph->src[kpart].n; ++iarr)
        if (!agraph->src[kpart].visited[iarr]) {
          inext = iarr;

          break;
        }

      size_t const jpart = agraph->src[kpart].itgt[inext];

      double xdiffij[BMM_NDIM];
      bmm_geom2d_diff(xdiffij, dem->part.x[jpart], dem->part.x[kpart]);

      double const gammai = $(bmm_swrap, double)(bmm_geom2d_dir(xdiffij), M_2PI);

      backdir = $(bmm_swrap, double)(gammai + M_PI, M_2PI);

      // fprintf(stderr, "Set back to %g (%zu).\n", backdir, inext);
    } else {
      double cand = M_2PI;
      for (size_t iarr = 0; iarr < agraph->src[kpart].n; ++iarr) {
        size_t jpart = agraph->src[kpart].itgt[iarr];

        double xdiffij[BMM_NDIM];
        bmm_geom2d_diff(xdiffij, dem->part.x[jpart], dem->part.x[kpart]);

        double const gammai = $(bmm_swrap, double)(bmm_geom2d_dir(xdiffij), M_2PI);

        if (!agraph->src[kpart].visited[iarr]) {
          // fprintf(stderr, "Unvisited angle from %zu to %zu of %g (%zu).\n", kpart, jpart, gammai, iarr);

          double const epsilonster = 1.0e-6;
          double const ncand = $(bmm_uwrap, double)(bmm_geom2d_dir(xdiffij) - backdir - epsilonster, M_2PI);
          if (ncand < cand) {
            cand = ncand;
            inext = iarr;
          }
        }
          // else fprintf(stderr, "Visited angle from %zu to %zu of %g (%zu).\n", kpart, jpart, gammai, iarr);
      }

      size_t const jpart = agraph->src[kpart].itgt[inext];

      double xdiffij[BMM_NDIM];
      bmm_geom2d_diff(xdiffij, dem->part.x[jpart], dem->part.x[kpart]);

      double const gammai = $(bmm_swrap, double)(bmm_geom2d_dir(xdiffij), M_2PI);

      backdir = $(bmm_swrap, double)(gammai + M_PI, M_2PI);

      // fprintf(stderr, "Set back to %g (%zu).\n", backdir, inext);
    }

    // fprintf(stderr, "Moved %zu -> %zu.\n", kpart, agraph->src[kpart].itgt[inext]);
    size_t const jpart = agraph->src[kpart].itgt[inext];

    agraph->src[kpart].visited[inext] = true;

    kpart = jpart;

    goto morer;
  }

  // Export.
  for (size_t iface = 0; iface < nface; ++iface) {
    for (size_t ivert = 0; ivert < faces->poly[iface].n; ++ivert) {
      size_t const ipart = faces->poly[iface].ivert[ivert];

      if (fprintf(stream, "%s%zu", ivert == 0 ? "" : " ", ipart) < 0) {
        BMM_TLE_STDS();

        return false;
      }
    }

    if (fprintf(stream, "\n") < 0) {
      BMM_TLE_STDS();

      return false;
    }
  }

  free(agraph);
  free(faces);

  if (fclose(stream) != 0) {
    BMM_TLE_STDS();

    return false;
  }

  return true;
}

// TODO These procedures semisuck, just like everything else here,
// but I thought I should still let you know this.

static double wkde_eval(double const *restrict const xarr,
    double const *restrict const warr,
    size_t const nsample, double const x, double const bandwidth) {
  double y = 0.0;

  for (size_t i = 0; i < nsample; i++)
    y += warr[i] * bmm_kernel_epan((x - xarr[i]) / bandwidth) / bandwidth;

  return y;
}

static void wkde_sample(double *restrict const yarr,
    double *restrict const yyarr,
    double const *restrict const xarr,
    double const *restrict const warr, size_t const nsample,
    double const bandwidth,
    size_t const narr, double const min, double const max) {
  double const step = (max - min) / (double) (narr - 1);

  for (size_t i = 0; i < narr; i++) {
    double const x = min + (double) i * step;

    yarr[i] = x;
    yyarr[i] = wkde_eval(xarr, warr, nsample, x, bandwidth);
  }
}

static double wkde_eval_sorted(double const *restrict const xarr,
    double const *restrict const warr,
    size_t const nsample, double const x, double const bandwidth,
    size_t const iwin, size_t const jwin) {
  double y = 0.0;

  for (size_t i = iwin; i < jwin; i++)
    y += warr[i] * bmm_kernel_epan((x - xarr[i]) / bandwidth) / bandwidth;

  return y;
}

static void wkde_sample_sorted(double *restrict const yarr,
    double *restrict const yyarr,
    double const *restrict const xarr,
    double const *restrict const warr, size_t const nsample,
    double const bandwidth,
    size_t const narr, double const min, double const max) {
  double const step = (max - min) / (double) (narr - 1);

  double const win = bandwidth;

  double x = min;

  // TODO Should use another `bsearch`-like function. This is also off by one.
  size_t iwin = 0;
  while (iwin < nsample && xarr[iwin] < x - win)
    ++iwin;
  size_t jwin = iwin;
  while (jwin < nsample && xarr[jwin] < x + win)
    ++jwin;

  for (size_t i = 0; i < narr; i++) {
    x = min + (double) i * step;

    yarr[i] = x;
    yyarr[i] = wkde_eval_sorted(xarr, warr, nsample, x, bandwidth, iwin, jwin);

    while (iwin < nsample && xarr[iwin] < x - win)
      ++iwin;
    while (jwin < nsample && xarr[jwin] < x + win)
      ++jwin;
  }
}

__attribute__ ((__nonnull__, __pure__))
static int compar(size_t const i, size_t const j, void *const cls) {
  double const *const *const rw = cls;

  return $(bmm_cmp, double)(rw[0][i], rw[0][j]);
}

__attribute__ ((__nonnull__))
static void swap(size_t const i, size_t const j, void *const cls) {
  double *const *const rw = cls;

  for (size_t k = 0; k < 2; ++k) {
    double const tmp = rw[k][i];
    rw[k][i] = rw[k][j];
    rw[k][j] = tmp;
  }
}

/// The call `bmm_dem_est_raddist(pr, pg, dem, nr)`
/// sets `pr` and `pg` of length `nr`
/// to the radial distribution function of the particles
/// in the simulation `dem`.
/// The simulation cell must be full for this to produce an accurate result.
__attribute__ ((__nonnull__))
bool bmm_dem_est_raddist(double *const pr, double *const pg,
    size_t const nbin, double const rmax,
    struct bmm_dem const *const dem) {
  // TODO This is bad and stupid.

  // double g[BMM_TRI(BMM_MPART)];
  size_t nmemb = $(bmm_tri, size_t)(dem->part.n);
  double *const w = malloc(nmemb * sizeof *w);
  if (w == NULL)
    return false;
  double *const r = malloc(nmemb * sizeof *r);
  if (r == NULL) {
    free(w);
    return false;
  }

// #define FUCK_WEIGHTS

  size_t i = 0;
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    if (bmm_dem_inside(dem, ipart))
      for (size_t jpart = ipart + 1; jpart < dem->part.n; ++jpart)
        if (bmm_dem_inside(dem, jpart)) {

          double const d = bmm_geom2d_cpdist(dem->part.x[ipart],
              dem->part.x[jpart], dem->opts.box.x, dem->opts.box.per);

          double const v0 = bmm_geom_ballsurf(d, 2);
          double const v = bmm_geom2d_shellvol(dem->part.x[ipart], d,
              dem->opts.box.x, dem->opts.box.per);

          r[i] = d;
#ifdef FUCK_WEIGHTS
          w[i] = 1.0;
#else
          w[i] = d == 0.0 ? 0.0 : v0 / v;
#endif
          ++i;
        }

  nmemb = i;

  // We convert `w` from "importance weights"
  // to "frequency weights" or "analytic weights".
  double wsum = $(bmm_sum, double)(w, nmemb);
  for (size_t i = 0; i < nmemb; ++i)
    w[i] /= wsum;

  double const bw = bmm_ival_midpoint(dem->opts.part.rnew) / 8.0;

  double *rw[] = {r, w};
  bmm_hsort_cls(nmemb, compar, swap, rw);
  wkde_sample_sorted(pr, pg, r, w, nmemb, bw, nbin, 0.0, rmax);
  // wkde_sample(pr, pg, r, w, nmemb, bw, nbin, 0.0, rmax);

  double const dr = rmax / (double) nbin;

  // Ideal gas.
  // for (size_t i = 0; i < nbin; ++i)
  //   pg[i] = bmm_geom_ballsurf(pr[i], BMM_NDIM);
  double total = 0.0;
  for (size_t i = 0; i < nbin; ++i)
    total += pg[i] * dr;
  for (size_t i = 0; i < nbin; ++i)
    pg[i] /= total;
  // It is a proper pdf now.
  for (size_t i = 0; i < nbin; ++i)
    pg[i] = pr[i] == 0.0 ? 0.0 : pg[i] *
      (bmm_geom_ballvol(rmax, BMM_NDIM) / bmm_geom_ballsurf(pr[i], BMM_NDIM));

  free(r);
  free(w);

  return true;
}

// TODO Complete energy estimators.

__attribute__ ((__nonnull__, __pure__))
double bmm_dem_est_epotext(struct bmm_dem const *const dem) {
  double e = 0.0;

  switch (dem->ext.tag) {
    case BMM_DEM_EXT_HARM:
      for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
        e += dem->ext.params.harm.k * $(bmm_power, double)
          (dem->opts.box.x[1] / 2.0 - dem->part.x[ipart][1], 2);

      break;
    case BMM_DEM_EXT_GRAVY:
      break;
    case BMM_DEM_EXT_DRIVE:
      break;
  }

  return (1.0 / 2.0) * e;
}

__attribute__ ((__nonnull__, __pure__))
double bmm_dem_est_eklin(struct bmm_dem const *const dem) {
  double e = 0.0;

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    for (size_t idim = 0; idim < BMM_NDIM; ++idim)
      e += dem->part.m[ipart] * $(bmm_power, double)(dem->part.v[ipart][idim], 2);

  return (1.0 / 2.0) * e;
}

__attribute__ ((__nonnull__, __pure__))
double bmm_dem_est_ekrot(struct bmm_dem const *const dem) {
  double e = 0.0;

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    e += dem->cache.j[ipart] * $(bmm_power, double)(dem->part.omega[ipart], 2);

  return (1.0 / 2.0) * e;
}

/// The call `bmm_dem_est_mass(dem)`
/// returns the total mass of the particles
/// in the simulation `dem`.
__attribute__ ((__nonnull__, __pure__))
double bmm_dem_est_mass(struct bmm_dem *const dem) {
  double m = 0.0;

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    m += dem->part.m[ipart];

  return m;
}

/// The call `bmm_dem_est_center(pxcenter, dem)`
/// sets `pxcenter` to the center of the bounding box
/// in the simulation `dem`.
__attribute__ ((__nonnull__))
void bmm_dem_est_center(double *const pxcenter,
    struct bmm_dem *const dem) {
  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    pxcenter[idim] = dem->opts.box.x[idim] / 2.0;
}

/// The call `bmm_dem_est_com(pxcom, dem)`
/// sets `pxcom` to the center of mass of the particles
/// in the simulation `dem`.
__attribute__ ((__nonnull__))
void bmm_dem_est_com(double *const pxcom,
    struct bmm_dem *const dem) {
  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    pxcom[idim] = 0.0;

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    for (size_t idim = 0; idim < BMM_NDIM; ++idim)
      pxcom[idim] += dem->part.m[ipart] * dem->part.x[ipart][idim];

  double const m = bmm_dem_est_mass(dem);

  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    pxcom[idim] /= m;
}

/// The call `bmm_dem_est_cor(dem)`
/// returns the mean coefficient of restitution of the particles
/// in the simulation `dem`.
/// The result only applies to the linear dashpot model and
/// even then it is a bit wrong.
__attribute__ ((__deprecated__, __nonnull__, __pure__))
double bmm_dem_est_cor(struct bmm_dem const *const dem) {
  double e = 0.0;

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    double const mred =
      dem->part.m[ipart] * dem->part.m[ipart] /
      (dem->part.m[ipart] + dem->part.m[ipart]);
    e += exp(-M_PI * dem->pair[BMM_DEM_CT_WEAK].norm.params.dashpot.gamma / (2.0 * mred) /
        sqrt(dem->pair[BMM_DEM_CT_WEAK].norm.params.dashpot.k / mred -
          $(bmm_power, double)(dem->pair[BMM_DEM_CT_WEAK].norm.params.dashpot.gamma / (2.0 * mred), 2)));
  }

  return e / (double) dem->part.n;
}

void bmm_dem_opts_def(struct bmm_dem_opts *const opts) {
  // This is here just to help Valgrind and cover up my mistakes.
  (void) memset(opts, 0, sizeof *opts);

  opts->verbose = false;

  opts->trap.enabled = false;
#ifdef _GNU_SOURCE
  opts->trap.mask = FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW;
#else
  opts->trap.mask = 0;
#endif

  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    opts->box.x[idim] = 1.0;

  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    opts->box.per[idim] = false;

  opts->time.istab = 1000;

  opts->part.y = 1.0;
  opts->part.nu = 0.5;
  opts->part.rnew[0] = 1.0;
  opts->part.rnew[1] = 1.0;

  opts->cont.ccrcont = 1.2;
  opts->cont.cshcont = 1.0;
  opts->cont.ktens = 1.0;
  opts->cont.dktens = 1.0;
  opts->cont.kshear = 1.0;
  opts->cont.dkshear = 1.0;

  opts->script.n = 0;

  opts->comm.dt = 1.0;
  opts->comm.flip = true;
  opts->comm.flop = true;
  opts->comm.flap = true;
  opts->comm.flup = true;

  for (size_t idim = 0; idim < BMM_NDIM; ++idim) {
    opts->cache.ncell[idim] = 5;

    // TODO Use this somewhere.
    /*
    dynamic_assert(opts->box.x[idim] / (double) (opts->cache.ncell[idim] * 2) >
        opts->part.rnew[1],
        "Neighbor cells too small");
    */
  }

  opts->cache.dcutoff = (double) INFINITY;
  for (size_t idim = 0; idim < BMM_NDIM; ++idim) {
    opts->cache.dcutoff = fmin(opts->cache.dcutoff,
        opts->box.x[idim] / (double) (opts->cache.ncell[idim] - 2));

    // TODO Use this somewhere.
    /*
    dynamic_assert(opts->cache.dcutoff <= opts->box.x[idim] /
        (double) ((opts->cache.ncell[idim] - 2) * 2),
        "Neighbor cells too small");
    */
  }
}

void bmm_dem_def(struct bmm_dem *const dem,
    struct bmm_dem_opts const *const opts) {
  // This is here just to help Valgrind and cover up my mistakes.
  (void) memset(dem, 0, sizeof *dem);

  dem->opts = *opts;

  dem->trap.remask = 0;

  dem->integ.tag = BMM_DEM_INTEG_EULER;
  dem->integ.tag = BMM_DEM_INTEG_TAYLOR;
  dem->integ.tag = BMM_DEM_INTEG_VELVET;
  dem->cache.tag = BMM_DEM_CACHE_NONE;
  dem->cache.tag = BMM_DEM_CACHE_NEIGH;
  dem->ext.tag = BMM_DEM_EXT_NONE;
  dem->amb.tag = BMM_DEM_AMB_FAXEN;
  dem->pair[BMM_DEM_CT_WEAK].norm.tag = BMM_DEM_NORM_KV;
  dem->pair[BMM_DEM_CT_WEAK].norm.tag = BMM_DEM_NORM_BSHP;
  dem->pair[BMM_DEM_CT_WEAK].tang.tag = BMM_DEM_TANG_HW;
  dem->pair[BMM_DEM_CT_WEAK].tang.tag = BMM_DEM_TANG_CS;
  dem->pair[BMM_DEM_CT_WEAK].torque.tag = BMM_DEM_TORQUE_SOFT;
  dem->pair[BMM_DEM_CT_WEAK].torque.tag = BMM_DEM_TORQUE_HARD;
  dem->yield.tag = BMM_DEM_YIELD_RANKINE;
  dem->yield.tag = BMM_DEM_YIELD_TRESCA;
  dem->cont_dep.tag = BMM_DEM_LINK_SPRING;
  dem->cont_dep.tag = BMM_DEM_LINK_BEAM;

  dem->amb.params.creeping.mu = 1.0e-3;

  switch (dem->pair[BMM_DEM_CT_WEAK].norm.tag) {
    case BMM_DEM_NORM_KV:
      dem->pair[BMM_DEM_CT_WEAK].norm.params.dashpot.k = 1.0e+6;
      dem->pair[BMM_DEM_CT_WEAK].norm.params.dashpot.gamma = 1.0e+3;

      break;
    case BMM_DEM_NORM_BSHP:
      dem->pair[BMM_DEM_CT_WEAK].norm.params.viscoel.a = 4.0e-2;

      break;
  }

  switch (dem->pair[BMM_DEM_CT_WEAK].tang.tag) {
    case BMM_DEM_TANG_HW:
      dem->pair[BMM_DEM_CT_WEAK].tang.params.hw.gamma = 1.0e+4;
      dem->pair[BMM_DEM_CT_WEAK].tang.params.hw.mu = 0.5;

      break;
    case BMM_DEM_TANG_CS:
      dem->pair[BMM_DEM_CT_WEAK].tang.params.cs.kappa = 1.0e+9;
      dem->pair[BMM_DEM_CT_WEAK].tang.params.cs.mu = 0.5;

      break;
  }

  switch (dem->cont_dep.tag) {
    case BMM_DEM_LINK_SPRING:

      break;
    case BMM_DEM_LINK_BEAM:

      break;
  }

  switch (dem->yield.tag) {
    case BMM_DEM_YIELD_RANKINE:
      dem->yield.params.rankine.sigmacrit = 1.0e+8;

      break;
    case BMM_DEM_YIELD_TRESCA:
      dem->yield.params.tresca.taucrit = 1.0e+8;

      break;
  }

  dem->time.t = 0.0;
  dem->time.istep = 0;

  dem->part.n = 0;
  dem->part.lnew = 0;

  for (size_t ipart = 0; ipart < nmembof(dem->pair[BMM_DEM_CT_STRONG].cont.src); ++ipart)
    dem->pair[BMM_DEM_CT_STRONG].cont.src[ipart].n = 0;

  dem->script.i = 0;
  dem->script.tprev = 0.0;

  dem->comm.tprev = 0.0;

  dem->cache.stale = false;
  dem->cache.i = 0;
  dem->cache.tpart = 0.0;
  dem->cache.tprev = 0.0;

  for (size_t icell = 0; icell < nmembof(dem->cache.part); ++icell)
    dem->cache.part[icell].n = 0;
}

// TODO Relocate these.

static bool msg_write(void const *buf, size_t const n,
    __attribute__ ((__unused__)) void *const ptr) {
  return bmm_io_writeout(buf, n);
}

size_t bmm_dem_sniff_size(struct bmm_dem const *const dem,
    enum bmm_msg_num const num) {
  switch (num) {
    case BMM_MSG_NUM_ISTEP:
      return sizeof dem->time;
    case BMM_MSG_NUM_OPTS:
      return sizeof dem->opts;
    case BMM_MSG_NUM_NEIGH:
      return sizeof dem->cache + sizeof dem->pair;
    case BMM_MSG_NUM_PARTS:
      return sizeof dem->part.n + sizeof dem->part;
    case BMM_MSG_NUM_EST:
      return sizeof dem->est;
  }

  dynamic_assert(false, "Unsupported message number");
}

bool bmm_dem_puts_stuff(struct bmm_dem const *const dem,
    enum bmm_msg_num const num) {
  switch (num) {
    case BMM_MSG_NUM_ISTEP:
      return msg_write(&dem->time, sizeof dem->time, NULL);
    case BMM_MSG_NUM_OPTS:
      return msg_write(&dem->opts, sizeof dem->opts, NULL);
    case BMM_MSG_NUM_NEIGH:
      return msg_write(&dem->cache, sizeof dem->cache, NULL) &&
        msg_write(&dem->pair, sizeof dem->pair, NULL);
    case BMM_MSG_NUM_PARTS:
      return msg_write(&dem->part.n, sizeof dem->part.n, NULL) &&
        msg_write(&dem->part, sizeof dem->part, NULL);
    case BMM_MSG_NUM_EST:
      return msg_write(&dem->est, sizeof dem->est, NULL);
  }

  dynamic_assert(false, "Unsupported message number");
}

bool bmm_dem_puts(struct bmm_dem const *const dem,
    enum bmm_msg_num const num) {
  struct bmm_msg_spec spec;
  bmm_msg_spec_def(&spec);
  spec.msg.size = bmm_dem_sniff_size(dem, num) + BMM_MSG_NUMSIZE;

  return bmm_msg_spec_write(&spec, msg_write, NULL) &&
    bmm_msg_num_write(&num, msg_write, NULL) &&
    bmm_dem_puts_stuff(dem, num);
}

bool bmm_dem_cache_expired(struct bmm_dem const *const dem) {
  // TODO Use `dem->opts.part.rnew[1]` instead of `dem->part.r[ipart]`.
  double const r = dem->opts.cache.dcutoff / 2.0;

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    if (bmm_geom2d_cpdist2(dem->part.x[ipart], dem->cache.x[ipart],
          dem->opts.box.x, dem->opts.box.per) >=
        $(bmm_power, double)(r - dem->part.r[ipart], 2))
      return true;

  return false;
}

__attribute__ ((__nonnull__))
static void bmm_dem_script_balance(struct bmm_dem *const dem) {
  double xcenter[BMM_NDIM];
  bmm_dem_est_center(xcenter, dem);

  double xcom[BMM_NDIM];
  bmm_dem_est_com(xcom, dem);

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    for (size_t idim = 0; idim < BMM_NDIM; ++idim) {
      dem->part.x[ipart][idim] += xcenter[idim] - xcom[idim];

      if (dem->opts.box.per[idim])
        dem->part.x[ipart][idim] = $(bmm_uwrap, double)(dem->part.x[ipart][idim],
            dem->opts.box.x[idim]);
    }
}

__attribute__ ((__nonnull__))
static bool bmm_dem_script_create_hc(struct bmm_dem *const dem) {
  double const etahc = bmm_geom_ballvol(0.5, BMM_NDIM);
  double const vhc = $(bmm_prod, double)(dem->opts.box.x, BMM_NDIM);
  double const eta = dem->opts.script.params[dem->script.i].create.eta;
  double const v = vhc * (etahc / eta);

  double const vlim = vhc * eta;
  // TODO Overlap factor.
  // double const rspace = bmm_ival_midpoint(dem->opts.part.rnew);
  double const rspace = dem->opts.part.rnew[1];

  double x[BMM_NDIM];
  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    x[idim] = 0.0;

  double vnow = 0.0;

  for ever {
    double const r = bmm_random_get(dem->rng, dem->opts.part.rnew);
    double const v = bmm_geom_ballvol(r, BMM_NDIM);

    double const vnext = vnow + v;

    if (vnext <= vlim) {
      size_t ipart = bmm_dem_addpart(dem);
      if (ipart == SIZE_MAX)
        return false;

      dem->part.r[ipart] = r;
      dem->part.m[ipart] = dem->opts.part.rho * v;

      for (size_t idim = 0; idim < BMM_NDIM; ++idim)
        dem->part.x[ipart][idim] = x[idim] + rspace;

      x[0] += 2.0 * rspace;

      if (x[0] + 2.0 * rspace >= dem->opts.box.x[0]) {
        x[0] = 0.0;
        x[1] += 2.0 * rspace;
      }

      vnow = vnext;
    } else
      break;
  }

  return true;
}

__attribute__ ((__nonnull__))
static bool bmm_dem_script_create_hex(struct bmm_dem *const dem) {
  double const etahc = bmm_geom_ballvol(0.5, BMM_NDIM);
  double const vhc = $(bmm_prod, double)(dem->opts.box.x, BMM_NDIM);
  double const eta = dem->opts.script.params[dem->script.i].create.eta;
  double const v = vhc * (etahc / eta);

  double const vlim = vhc * eta;
  // double const rspace = bmm_ival_midpoint(dem->opts.part.rnew);
  double const rspace = dem->opts.part.rnew[1];
  double const dspace = (sqrt(3.0) / 2.0) * rspace;

  double x[BMM_NDIM];
  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    x[idim] = 0.0;

  double vnow = 0.0;
  bool parity = false;

  for ever {
    double const r = bmm_random_get(dem->rng, dem->opts.part.rnew);
    double const v = bmm_geom_ballvol(r, BMM_NDIM);

    double const vnext = vnow + v;

    if (vnext <= vlim) {
      size_t ipart = bmm_dem_addpart(dem);
      if (ipart == SIZE_MAX)
        return false;

      dem->part.r[ipart] = r;
      dem->part.m[ipart] = dem->opts.part.rho * v;

      for (size_t idim = 0; idim < BMM_NDIM; ++idim)
        dem->part.x[ipart][idim] = x[idim] + rspace;

      x[0] += 2.0 * rspace;

      if (x[0] + 2.0 * rspace >= dem->opts.box.x[0]) {
        x[0] = parity ? 0.0 : rspace;
        x[1] += 2.0 * dspace;
        parity = !parity;
      }

      vnow = vnext;
    } else
      break;
  }

  return true;
}

// TODO Remove or polish this test diamond.
__attribute__ ((__nonnull__))
static bool bmm_dem_script_create_testpile(struct bmm_dem *const dem) {
  double const etahc = bmm_geom_ballvol(0.5, BMM_NDIM);
  double const vhc = $(bmm_prod, double)(dem->opts.box.x, BMM_NDIM);
  double const eta = dem->opts.script.params[dem->script.i].test.eta;
  double const v = vhc * (etahc / eta);

  double const vlim = vhc * eta;
  double const rspace = bmm_ival_midpoint(dem->opts.part.rnew);
  // double const rspace = dem->opts.part.rnew[1];
  double const dspace = (sqrt(3.0) / 2.0) * rspace;

  double x[BMM_NDIM];
  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    x[idim] = 0.0;

  double vnow = 0.0;
  bool parity = false;

  for ever {
    double const r = bmm_random_get(dem->rng, dem->opts.part.rnew);
    double const v = bmm_geom_ballvol(r, BMM_NDIM);

    double const vnext = vnow + v;

    if (vnext <= vlim) {
      size_t ipart = bmm_dem_addpart(dem);
      if (ipart == SIZE_MAX)
        return false;

      dem->part.r[ipart] = r;
      dem->part.m[ipart] = dem->opts.part.rho * v;

      for (size_t idim = 0; idim < BMM_NDIM; ++idim)
        dem->part.x[ipart][idim] = x[idim] + rspace;

      x[0] += 2.0 * rspace;

      if (x[0] + 2.0 * rspace >= dem->opts.box.x[0] + 2.0 * rspace) {
        x[0] = parity ? 0.0 : rspace;
        x[1] += 2.0 * dspace;
        parity = !parity;
      }

      if (fabs(dem->part.x[ipart][0] - dem->opts.box.x[0] / 2.0) +
          fabs(dem->part.x[ipart][1] - dem->opts.box.x[1] / 2.0) >
          dem->opts.box.x[0] / 1.9 ||
          dem->part.x[ipart][1] < dem->opts.box.x[1] / 2.0 - rspace)
        bmm_dem_rempart(dem, ipart);

      if (fabs(dem->part.x[ipart][1] - dem->opts.box.x[1] / 2.0) < rspace)
        dem->part.role[ipart] = BMM_DEM_ROLE_FIXED;

      vnow = vnext;
    } else
      break;
  }

  return true;
}

__attribute__ ((__nonnull__))
static bool bmm_dem_script_create_testplane(struct bmm_dem *const dem) {
  double const etahc = bmm_geom_ballvol(0.5, BMM_NDIM);
  double const vhc = $(bmm_prod, double)(dem->opts.box.x, BMM_NDIM);
  double const eta = dem->opts.script.params[dem->script.i].test.eta;
  double const v = vhc * (etahc / eta);

  double const vlim = vhc * eta;
  double const rspace = bmm_ival_midpoint(dem->opts.part.rnew);
  // double const rspace = dem->opts.part.rnew[1];
  double const dspace = (sqrt(3.0) / 2.0) * rspace;

  double x[BMM_NDIM];
  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    x[idim] = 0.0;

  double vnow = 0.0;
  bool parity = false;

  for ever {
    double const r = bmm_random_get(dem->rng, dem->opts.part.rnew);
    double const v = bmm_geom_ballvol(r, BMM_NDIM);

    double const vnext = vnow + v;

    if (vnext <= vlim) {
      size_t ipart = bmm_dem_addpart(dem);
      if (ipart == SIZE_MAX)
        return false;

      dem->part.r[ipart] = r;
      dem->part.m[ipart] = dem->opts.part.rho * v;

      for (size_t idim = 0; idim < BMM_NDIM; ++idim)
        dem->part.x[ipart][idim] = x[idim] + rspace;

      x[0] += 2.0 * rspace / 4.0;

      if (x[0] + 2.0 * rspace >= dem->opts.box.x[0] + 2.0 * rspace) {
        break;
      }

      dem->part.role[ipart] = BMM_DEM_ROLE_FIXED;

      vnow = vnext;
    } else
      break;
  }

  return true;
}

// TODO Remove or polish this test beam.
__attribute__ ((__nonnull__))
static bool bmm_dem_script_create_testbeam(struct bmm_dem *const dem) {
  double const etahc = bmm_geom_ballvol(0.5, BMM_NDIM);
  double const vhc = $(bmm_prod, double)(dem->opts.box.x, BMM_NDIM);
  double const eta = dem->opts.script.params[dem->script.i].test.eta;
  double const v = vhc * (etahc / eta);

  double const vlim = vhc * eta;
  // double const rspace = bmm_ival_midpoint(dem->opts.part.rnew);
  double const rspace = dem->opts.part.rnew[1];
  double const dspace = (sqrt(3.0) / 2.0) * rspace;

  double x[BMM_NDIM];
  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    x[idim] = 0.0;

  double vnow = 0.0;
  bool parity = false;

  for ever {
    double const r = bmm_random_get(dem->rng, dem->opts.part.rnew);
    double const v = bmm_geom_ballvol(r, BMM_NDIM);

    double const vnext = vnow + v;

    if (vnext <= vlim) {
      size_t ipart = bmm_dem_addpart(dem);
      if (ipart == SIZE_MAX)
        return false;

      dem->part.r[ipart] = r;
      dem->part.m[ipart] = dem->opts.part.rho * v;

      for (size_t idim = 0; idim < BMM_NDIM; ++idim)
        dem->part.x[ipart][idim] = x[idim] + rspace;

      x[0] += 2.0 * rspace;

      if (x[0] + 2.0 * rspace >= dem->opts.box.x[0]) {
        x[0] = parity ? 0.0 : rspace;
        x[1] += 2.0 * dspace;
        parity = !parity;
      }

      if (dem->part.x[ipart][1] > dem->opts.script.params[dem->script.i].test.layers * rspace ||
          dem->part.x[ipart][0] + 4.0 * rspace > dem->opts.script.params[dem->script.i].test.slices * rspace)
        bmm_dem_rempart(dem, ipart);

      if (dem->part.x[ipart][0] > dem->opts.box.x[0] - 4.0 * rspace)
        bmm_dem_rempart(dem, ipart);

      if (dem->part.x[ipart][0] < 4.0 * rspace)
        dem->part.role[ipart] = BMM_DEM_ROLE_FIXED;

      vnow = vnext;
    } else
      break;
  }

  for (size_t i = 0; i < 14; ++i) {
    size_t ipart = bmm_dem_addpart(dem);
    if (ipart == SIZE_MAX)
      return false;

    double const r = bmm_random_get(dem->rng, dem->opts.part.rnew);
    double const v = bmm_geom_ballvol(r, BMM_NDIM);

    dem->part.r[ipart] = r;
    dem->part.m[ipart] = dem->opts.part.rho * v;

    double a[2];
    a[0] = 0.0;
    a[1] = dem->opts.box.x[0];
    dem->part.x[ipart][0] = bmm_random_get(dem->rng, a);
    a[0] = 0.0;
    a[1] = sqrt(2.0 * dem->opts.box.x[1]);
    dem->part.x[ipart][1] =
      bmm_random_get(dem->rng, a) * bmm_random_get(dem->rng, a) -
      $(bmm_power, double)(bmm_ival_length(a), 2);

    dem->part.role[ipart] = BMM_DEM_ROLE_FIXED;
  }

  return true;
}

// TODO Bad!
__attribute__ ((__nonnull__))
static void bmm_dem_script_perturb(struct bmm_dem *const dem) {
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    for (size_t idim = 0; idim < BMM_NDIM; ++idim)
      dem->part.x[ipart][idim] += bmm_random_get(dem->rng,
          dem->opts.part.rnew) / 4.0;
}

__attribute__ ((__nonnull__))
static bool bmm_dem_script_create_gas(struct bmm_dem *const dem) {
  for (size_t ipart = 0; ipart < 64; ++ipart) {
    size_t const jpart = bmm_dem_addpart(dem);
    if (jpart == SIZE_MAX)
      return false;

    dem->part.r[jpart] = 0.03;
    dem->part.m[jpart] = 1.0;

    for (size_t idim = 0; idim < BMM_NDIM; ++idim)
      dem->part.x[jpart][idim] += gsl_rng_uniform(dem->rng) * dem->opts.box.x[idim];

      dem->part.omega[jpart] += (double) (rand() % 512 - 256);
  }

  return true;
}

__attribute__ ((__nonnull__))
static bool bmm_dem_script_create_couple(struct bmm_dem *const dem) {
  size_t jpart;

  jpart = bmm_dem_addpart(dem);
  if (jpart == SIZE_MAX)
    return false;

  double const scale = dem->opts.box.x[0];

  dem->part.r[jpart] = 0.03 * scale;
  dem->part.m[jpart] = dem->opts.part.rho * bmm_geom_ballvol(dem->part.r[jpart], BMM_NDIM);
  dem->part.x[jpart][0] += 0.45 * scale;
  dem->part.x[jpart][1] += 0.25 * scale;
  dem->part.v[jpart][0] += 3.0e+1 * scale;
  dem->part.omega[jpart] += 3.0e+4;

  jpart = bmm_dem_addpart(dem);
  if (jpart == SIZE_MAX)
    return false;

  dem->part.r[jpart] = 0.03 * scale;
  dem->part.m[jpart] = dem->opts.part.rho * bmm_geom_ballvol(dem->part.r[jpart], BMM_NDIM);
  dem->part.x[jpart][0] += 0.55 * scale;
  dem->part.x[jpart][1] += 0.25 * scale;
  dem->part.v[jpart][0] -= 3.0e+1 * scale;
  dem->part.omega[jpart] += 3.0e+4;

  return true;
}

__attribute__ ((__nonnull__))
static bool bmm_dem_script_create_triplet(struct bmm_dem *const dem) {
  size_t jpart;

  jpart = bmm_dem_addpart(dem);
  if (jpart == SIZE_MAX)
    return false;

  double const scale = dem->opts.box.x[0];

  dem->part.r[jpart] = 0.03 * scale;
  dem->part.m[jpart] = dem->opts.part.rho * bmm_geom_ballvol(dem->part.r[jpart], BMM_NDIM);
  dem->part.x[jpart][0] += 0.45 * scale;
  dem->part.x[jpart][1] += 0.25 * scale;
  dem->part.v[jpart][0] += 3.0e+1 * scale;
  dem->part.omega[jpart] += 3.0e+4;

  jpart = bmm_dem_addpart(dem);
  if (jpart == SIZE_MAX)
    return false;

  dem->part.r[jpart] = 0.03 * scale;
  dem->part.m[jpart] = dem->opts.part.rho * bmm_geom_ballvol(dem->part.r[jpart], BMM_NDIM);
  dem->part.x[jpart][0] += 0.55 * scale;
  dem->part.x[jpart][1] += 0.25 * scale;
  dem->part.v[jpart][0] -= 3.0e+1 * scale;
  dem->part.omega[jpart] += 3.0e+4;

  jpart = bmm_dem_addpart(dem);
  if (jpart == SIZE_MAX)
    return false;

  dem->part.r[jpart] = 0.02 * scale;
  dem->part.m[jpart] = dem->opts.part.rho * bmm_geom_ballvol(dem->part.r[jpart], BMM_NDIM);
  dem->part.x[jpart][0] += 0.5 * scale;
  dem->part.x[jpart][1] += 0.25 * scale;
  dem->part.v[jpart][0] -= 0.0e+1 * scale;
  dem->part.omega[jpart] += -3.0e+4;

  return true;
}

/// The call `bmm_dem_step(dem)`
/// advances the simulation `dem` by one step.
/// Make sure the simulation has not ended prior to the call
/// by calling `bmm_dem_script_ongoing` or `bmm_dem_script_trans`.
bool bmm_dem_step(struct bmm_dem *const dem) {
  switch (dem->opts.script.mode[dem->script.i]) {
    case BMM_DEM_MODE_CREATE_HC:
      if (!bmm_dem_script_create_hc(dem))
        return false;

      bmm_dem_script_perturb(dem);

      bmm_dem_script_balance(dem);

      break;
    case BMM_DEM_MODE_CREATE_HEX:
      if (!bmm_dem_script_create_hex(dem))
        return false;

      bmm_dem_script_perturb(dem);

      bmm_dem_script_balance(dem);

      break;
    case BMM_DEM_MODE_CREATE_PLANE:
      if (!bmm_dem_script_create_testplane(dem))
        return false;

      bmm_dem_script_perturb(dem);

      break;
    case BMM_DEM_MODE_CREATE_PILE:
      if (!bmm_dem_script_create_testpile(dem))
        return false;

      bmm_dem_script_perturb(dem);

      bmm_dem_script_balance(dem);

      break;
    case BMM_DEM_MODE_CREATE_BEAM:
      if (!bmm_dem_script_create_testbeam(dem))
        return false;

      bmm_dem_script_perturb(dem);

      bmm_dem_script_balance(dem);

      break;
    case BMM_DEM_MODE_CREATE_GAS:
      if (!bmm_dem_script_create_gas(dem))
        return false;

      bmm_dem_script_balance(dem);

      break;
    case BMM_DEM_MODE_CREATE_COUPLE:
      if (!bmm_dem_script_create_couple(dem))
        return false;

      bmm_dem_script_balance(dem);

      break;
    case BMM_DEM_MODE_CREATE_TRIPLET:
      if (!bmm_dem_script_create_triplet(dem))
        return false;

      bmm_dem_script_balance(dem);

      break;
    case BMM_DEM_MODE_GRAVY:
      dem->ext.tag = BMM_DEM_EXT_GRAVY;
      dem->ext.params.gravy.g =
        -dem->opts.script.params[dem->script.i].gravy.g;

      break;
    case BMM_DEM_MODE_SEDIMENT:
      dem->ext.tag = BMM_DEM_EXT_HARM;
      dem->ext.params.harm.k =
        dem->opts.script.params[dem->script.i].sediment.kcohes;

      break;
    case BMM_DEM_MODE_PRESET0:
      dem->opts.comm.dt *= 1.0e+1;
      dem->amb.tag = BMM_DEM_AMB_FAXEN;
      dem->pair[BMM_DEM_CT_WEAK].norm.tag = BMM_DEM_NORM_BSHP;
      dem->pair[BMM_DEM_CT_WEAK].norm.params.viscoel.a = 2.0e-2;
      dem->pair[BMM_DEM_CT_WEAK].tang.tag = BMM_DEM_TANG_NONE;
      dem->cont_dep.tag = BMM_DEM_LINK_SPRING;

      break;
    case BMM_DEM_MODE_CLIP:
      bmm_dem_script_clip(dem);

      break;
    case BMM_DEM_MODE_PRESET1:
      dem->opts.comm.dt /= 1.0e+1;
      dem->amb.tag = BMM_DEM_AMB_NONE;
      dem->pair[BMM_DEM_CT_WEAK].norm.tag = BMM_DEM_NORM_BSHP;
      dem->pair[BMM_DEM_CT_WEAK].tang.tag = BMM_DEM_TANG_HW;
      dem->pair[BMM_DEM_CT_WEAK].tang.tag = BMM_DEM_TANG_CS;
      dem->cont_dep.tag = BMM_DEM_LINK_BEAM;
      // dem->opts.part.y = 52.0e+9;

      break;
    case BMM_DEM_MODE_LINK:
      dem->ext.tag = BMM_DEM_EXT_NONE;

      if (!bmm_dem_link(dem))
        return false;

      // TODO Yoink test.
      // dem->part.omega[0] = 1.0e+6;

      break;
    case BMM_DEM_MODE_FAULT:
      bmm_dem_fault(dem);

      break;
    case BMM_DEM_MODE_SEPARATE:

      break;
    case BMM_DEM_MODE_GLUE:
      for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
        // TODO Parametrize the glue layer thickness.
        if (dem->part.x[ipart][1] < dem->opts.box.x[1] * 0.1)
          dem->part.role[ipart] = BMM_DEM_ROLE_FIXED;

      break;
    case BMM_DEM_MODE_CRUNCH:
      dem->ext.tag = BMM_DEM_EXT_DRIVE;

      static bool pristine = true;
      if (pristine) {
        pristine = false;

        for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
          if (dem->part.x[ipart][1] > dem->opts.box.x[1] * 0.9)
            dem->part.role[ipart] = BMM_DEM_ROLE_DRIVEN;

        for (size_t idim = 0; idim < BMM_NDIM; ++idim)
          dem->script.state.crunch.f[idim] = 0.0;
      }

      break;
    case BMM_DEM_MODE_MEASURE:

      break;
  }

  if (dem->cache.stale || bmm_dem_cache_expired(dem)) {
    if (!bmm_dem_cache_build(dem))
      return false;

    dem->cache.tprev = dem->time.t;
  }

  bmm_dem_predict(dem);
  bmm_dem_analyze(dem);
  bmm_dem_force(dem);
  bmm_dem_accel(dem);
  bmm_dem_correct(dem);

  if (dem->time.istep % dem->opts.time.istab == 0)
    bmm_dem_stab(dem);

  dem->time.t += dem->opts.script.dt[dem->script.i];
  ++dem->time.istep;

  return true;
}

static FILE *stream;

static bool pregarbage(struct bmm_dem const *const dem) {
  stream = fopen("garbage.data", "w");
  if (stream == NULL) {
    BMM_TLE_STDS();

    return false;
  }

  return true;
}

static bool garbage(struct bmm_dem const *const dem) {
  if (fprintf(stream, "%g %g %g %g\n",
        dem->time.t,
        dem->est.epotext,
        dem->est.mueff,
        dem->est.vdriv[0]) < 0) {
    BMM_TLE_STDS();

    return false;
  }

  return true;
}

static bool postgarbage(struct bmm_dem const *const dem) {
  if (fclose(stream) != 0) {
    BMM_TLE_STDS();

    return false;
  }

  return true;
}

// TODO This looks just like `bmm_dem_script_trans`.
bool bmm_dem_comm(struct bmm_dem *const dem) {
  // TODO Make a mechanism to automate retransmission of differences only.

  double const toff = dem->time.t - dem->comm.tprev - dem->opts.comm.dt;

  if (toff >= 0.0) {
    dem->comm.tprev = dem->time.t;

    // TODO Nope.
    static bool first = true;
    if (first) {
      if (!bmm_dem_puts(dem, BMM_MSG_NUM_OPTS))
        return false;

      first = false;
    }

    if (!bmm_dem_puts(dem, BMM_MSG_NUM_ISTEP))
      return false;

    if (!bmm_dem_puts(dem, BMM_MSG_NUM_NEIGH))
      return false;

    if (!bmm_dem_puts(dem, BMM_MSG_NUM_PARTS))
      return false;

    dem->est.epotext = bmm_dem_est_epotext(dem);
    dem->est.eklin = bmm_dem_est_eklin(dem);
    dem->est.ekrot = bmm_dem_est_ekrot(dem);

    if (!bmm_dem_puts(dem, BMM_MSG_NUM_EST))
      return false;

    if (!garbage(dem))
      BMM_TLE_EXTS(BMM_TLE_NUM_UNKNOWN, "Fucked up");

    // if (!extra_crap(dem))
    //   return false;
  }

  return true;
}

static bool rubbish(struct bmm_dem const *const dem) {
  FILE *const stream = fopen("rubbish.data", "w");
  if (stream == NULL) {
    BMM_TLE_STDS();

    return false;
  }

  size_t nbin = 512;
  double *const r = malloc(nbin * sizeof *r);
  double *const g = malloc(nbin * sizeof *g);
  dynamic_assert(r != NULL && g != NULL, "Allocated");

  double const rmax = bmm_fp_min(dem->opts.box.x, BMM_NDIM) / 2.0;

  if (!bmm_dem_est_raddist(r, g, nbin, rmax, dem))
    return false;

  for (size_t ibin = 0; ibin < nbin; ++ibin)
    if (fprintf(stream, "%g %g\n", r[ibin], g[ibin]) < 0) {
      BMM_TLE_STDS();

      return false;
    }

  // TODO Cut this.

  // Snip.

  double const v = $(bmm_prod, double)(dem->opts.box.x, BMM_NDIM);
  double const rho = (double) dem->part.n / v;
  double const dr = rmax / (double) nbin;
  double a = 0.0;

  for (size_t ibin = 0; ibin < nbin; ++ibin)
    a += (g[ibin] * log(g[ibin] + 1.0e-9) - (g[ibin] - 1.0)) * dr;

  double const sk = -(rho / 2.0) * a;

  if (fprintf(stderr, "s / k = %g\n", sk) < 0) {
    BMM_TLE_STDS();

    return false;
  }

  // Snap.

  free(g);
  free(r);

  if (fclose(stream) != 0) {
    BMM_TLE_STDS();

    return false;
  }

  return true;
}

static double abserr(size_t const i, double const z, void const *const cls) {
  double const *const t = cls;

  return z + fabs(t[i]);
}

bool bmm_dem_report(struct bmm_dem const *const dem) {
  if (dem->opts.verbose) {
    if (fprintf(stderr, "Time Error: %g\n",
          $(bmm_foldl_cls, double)(dem->opts.script.n,
            abserr, 0.0, dem->script.toff)) < 0)
      return false;
  }

  return true;
}

static bool bmm_dem_run_(struct bmm_dem *const dem) {
  int const sigs[] = {SIGINT, SIGQUIT, SIGTERM, SIGPIPE};
  if (bmm_sig_register(sigs, nmembof(sigs)) != SIZE_MAX) {
    BMM_TLE_STDS();

    return false;
  }

  if (!pregarbage(dem))
    BMM_TLE_EXTS(BMM_TLE_NUM_UNKNOWN, "Fucked up");

  for ever {
    int signum;
    if (bmm_sig_use(&signum))
      switch (signum) {
        case SIGINT:
        case SIGQUIT:
        case SIGTERM:
        case SIGPIPE:
          BMM_TLE_EXTS(BMM_TLE_NUM_ASYNC, "Interrupted");

          return false;
      }

    if (!bmm_dem_script_ongoing(dem))
      return true;

    if (!bmm_dem_comm(dem))
      return false;

    if (!bmm_dem_step(dem))
      return false;

    if (!bmm_dem_script_trans(dem))
      return true;
  }

  return true;
}

bool bmm_dem_trap_on(struct bmm_dem *const dem) {
  if (dem->opts.trap.enabled) {
#ifdef _GNU_SOURCE
    dem->trap.remask = feenableexcept(dem->opts.trap.mask);
    if (dem->trap.remask == -1) {
      BMM_TLE_STDS();

      return false;
    }
#else
    BMM_TLE_EXTS(BMM_TLE_NUM_UNSUPP, "Trapping exceptions unsupported");

    return false;
#endif
  }

  return true;
}

bool bmm_dem_trap_off(struct bmm_dem *const dem) {
  if (dem->opts.trap.enabled) {
#ifdef _GNU_SOURCE
    if (feenableexcept(dem->trap.remask) == -1) {
      BMM_TLE_STDS();

      return false;
    }
#else
    BMM_TLE_EXTS(BMM_TLE_NUM_UNSUPP, "Trapping exceptions unsupported");

    return false;
#endif
  }

  return true;
}

bool bmm_dem_run(struct bmm_dem *const dem) {
// #define POST_DEBUG
#ifndef POST_DEBUG
  bmm_dem_trap_on(dem);

  bool const run = bmm_dem_run_(dem);
  bool const report = bmm_dem_report(dem);

  FILE *const stream = fopen("a.out", "w");
  if (stream == NULL) {
    BMM_TLE_STDS();
    abort();
  }
  if (fwrite(dem, sizeof *dem, 1, stream) != 1) {
    BMM_TLE_STDS();
    abort();
  }
  if (fclose(stream) != 0) {
    BMM_TLE_STDS();
    abort();
  }

  bmm_dem_trap_off(dem);

  if (!postgarbage(dem))
    BMM_TLE_EXTS(BMM_TLE_NUM_UNKNOWN, "Fucked up");

  if (!rubbish(dem))
    BMM_TLE_EXTS(BMM_TLE_NUM_UNKNOWN, "Fucked up");
#else
  bool const run = true;
  bool const report = true;

  // Uh oh!
  gsl_rng *const rng = dem->rng;
  FILE *const strim = fopen("a.out", "r");
  if (strim == NULL) {
    BMM_TLE_STDS();
    abort();
  }
  if (fread(dem, sizeof *dem, 1, strim) != 1) {
    BMM_TLE_STDS();
    abort();
  }
  if (fclose(strim) != 0) {
    BMM_TLE_STDS();
    abort();
  }
  dem->rng = rng;
#endif

  if (!export_x(dem) || !export_r(dem) || !export_c(dem) ||
      !export_f(dem) || !export_p(dem))
    BMM_TLE_EXTS(BMM_TLE_NUM_UNKNOWN, "Fucked up big");

  return run && report;
}

static bool bmm_dem_run_with_(struct bmm_dem *const dem) {
  gsl_rng_type const *const t = gsl_rng_env_setup();
  if (t == NULL) {
    BMM_TLE_STDS();

    return false;
  }

  dem->rng = gsl_rng_alloc(t);
  if (dem->rng == NULL) {
    BMM_TLE_STDS();

    return false;
  }

  bool const result = bmm_dem_run(dem);

  gsl_rng_free(dem->rng);

  return result;
}

bool bmm_dem_run_with(struct bmm_dem_opts const *const opts) {
  struct bmm_dem *const dem = malloc(sizeof *dem);
  if (dem == NULL) {
    BMM_TLE_STDS();

    return false;
  }

  bmm_dem_def(dem, opts);

  bool const result = bmm_dem_run_with_(dem);

  free(dem);

  return result;
}
