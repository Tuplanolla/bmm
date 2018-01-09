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

// Apologies for the horrible mess that this file became.

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

void bmm_dem_opts_set_rnew(struct bmm_dem_opts *const opts,
    double const *const rnew) {
  double const leeway = 4.0;

  opts->part.rnew[0] = rnew[0];
  opts->part.rnew[1] = rnew[1];

  opts->cache.dcutoff = leeway * rnew[1];

  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    opts->cache.ncell[idim] = $(bmm_max, size_t)(3,
      (size_t) (opts->box.x[idim] / opts->cache.dcutoff)) +
      (opts->box.per[idim] ? 0 : 2);

  for (size_t idim = 0; idim < BMM_NDIM; ++idim) {
    dynamic_assert(opts->box.x[idim] / (double) opts->cache.ncell[idim] >
        2.0 * opts->part.rnew[1],
        "Neighbor cells too small");

    dynamic_assert(opts->cache.dcutoff <= opts->box.x[idim] /
        (double) ((opts->cache.ncell[idim] - (opts->box.per[idim] ? 0 : 2))),
        "Neighbor cells too small");
  }
}

__attribute__ ((__nonnull__))
void bmm_dem_est_vdc(double *const pv, struct bmm_dem const *const dem) {
  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    pv[idim] = 0.0;

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    if (dem->part.role[ipart] == BMM_DEM_ROLE_DRIVEN)
      for (size_t idim = 0; idim < BMM_NDIM; ++idim)
        pv[idim] += dem->part.v[ipart][idim];

  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    pv[idim] /= (double) dem->script.state.crunch.ndrive;
}

__attribute__ ((__nonnull__, __pure__))
double bmm_dem_est_epotext(struct bmm_dem const *const dem) {
  double e = 0.0;

  switch (dem->ext.tag) {
    case BMM_DEM_EXT_HARM:
      for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
        e += (1.0 / 2.0) * dem->ext.params.harm.k * $(bmm_power, double)
          (dem->opts.box.x[1] / 2.0 - dem->part.x[ipart][1], 2);

      break;
    case BMM_DEM_EXT_GRAVY:
      for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
        e += dem->ext.params.gravy.g * dem->part.m[ipart] *
          (dem->opts.box.x[1] / 2.0 - dem->part.x[ipart][1]);

      break;
    case BMM_DEM_EXT_DRIVE:
      break;
  }

  return e;
}

__attribute__ ((__nonnull__, __pure__))
double bmm_dem_est_eklin(struct bmm_dem const *const dem) {
  double e = 0.0;

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    for (size_t idim = 0; idim < BMM_NDIM; ++idim)
      e += dem->part.m[ipart] *
        $(bmm_power, double)(dem->part.v[ipart][idim], 2);

  return (1.0 / 2.0) * e;
}

__attribute__ ((__nonnull__, __pure__))
double bmm_dem_est_ekrot(struct bmm_dem const *const dem) {
  double e = 0.0;

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    e += dem->cache.j[ipart] * $(bmm_power, double)(dem->part.omega[ipart], 2);

  return (1.0 / 2.0) * e;
}

__attribute__ ((__nonnull__, __pure__))
double bmm_dem_est_econt_one(struct bmm_dem const *const dem,
    enum bmm_dem_ct const ict, size_t const ipart, size_t const icont, size_t const jpart) {
  double e = 0.0;

  double xdiffij[BMM_NDIM];
  bmm_geom2d_cpdiff(xdiffij, dem->part.x[jpart], dem->part.x[ipart],
      dem->opts.box.x, dem->opts.box.per);

  double const d2 = bmm_geom2d_norm2(xdiffij);
  if (d2 == 0.0)
    return e;

  double const ri = dem->part.r[ipart];
  double const rj = dem->part.r[jpart];
  double const r = ri + rj;
  double const r2 = $(bmm_power, double)(r, 2);
  if (!dem->pair[ict].cohesive)
    if (d2 > r2)
      return e;

  double const d = sqrt(d2);

  double xnormij[BMM_NDIM];
  bmm_geom2d_scale(xnormij, xdiffij, 1.0 / d);

  double xtangij[BMM_NDIM];
  bmm_geom2d_rperp(xtangij, xnormij);

  double vdiffij[BMM_NDIM];
  bmm_geom2d_diff(vdiffij, dem->part.v[jpart], dem->part.v[ipart]);

  double const xi = r - d;
  double const vnormij = -bmm_geom2d_dot(vdiffij, xnormij);
  double const reff = $(bmm_resum2, double)(ri, rj);
  double const dt = dem->opts.script.dt[dem->script.i];

  double fnorm = 0.0;

  switch (dem->pair[ict].norm.tag) {
    case BMM_DEM_NORM_KV:
      if (dem->pair[ict].cohesive) {
        // double const r = ri + rj;
        // Maybe.
        double const r = dem->pair[ict].cont.src[ipart].drest[icont];
        double const xi = r - d;

        fnorm = dem->pair[ict].norm.params.dashpot.k * xi +
          dem->pair[ict].norm.params.dashpot.gamma * vnormij;

        e += (1.0 / 2.0) * dem->pair[ict].norm.params.dashpot.k * $(bmm_power, double)(xi, 2);
      } else {
        fnorm = $(bmm_max, double)(0.0,
            dem->pair[ict].norm.params.dashpot.k * xi +
            dem->pair[ict].norm.params.dashpot.gamma * vnormij);

        if (0.0 > dem->pair[ict].norm.params.dashpot.k * xi +
            dem->pair[ict].norm.params.dashpot.gamma * vnormij) {
        } else {
          e += (1.0 / 2.0) * dem->pair[ict].norm.params.dashpot.k * $(bmm_power, double)(xi, 2);
        }
      }

      break;
    case BMM_DEM_NORM_BSHP:
      {
        double const more = dem->est.bshpp * sqrt(reff * xi);

        if (dem->pair[ict].cohesive) {
          fnorm = more * (xi + dem->pair[ict].norm.params.viscoel.a * vnormij);

          e += (2.0 / 5.0) * more * $(bmm_power, double)(xi, 2);
        } else {
          fnorm = $(bmm_max, double)(0.0,
              more * (xi + dem->pair[ict].norm.params.viscoel.a * vnormij));

          if (0.0 > more * (xi + dem->pair[ict].norm.params.viscoel.a * vnormij)) {
          } else {
            e += (2.0 / 5.0) * more * $(bmm_power, double)(xi, 2);
          }
        }
      }

      break;
  }

  switch (dem->pair[ict].tang.tag) {
    case BMM_DEM_TANG_HW:

      break;
    case BMM_DEM_TANG_CS:
      {
        double const lambdaij = bmm_geom2d_dir(xdiffij);
        double const lambdaji = $(bmm_swrap, double)(lambdaij + M_PI, M_2PI);
        double const phii = dem->part.phi[ipart];
        double const phij = dem->part.phi[jpart];
        double const psii = phii - lambdaij;
        double const psij = phij - lambdaji;

        double const dpsii = $(bmm_swrap, double)(psii - dem->pair[ict].cont.src[ipart].psirest[icont][BMM_DEM_END_TAIL], M_2PI);
        double const dpsij = $(bmm_swrap, double)(psij - dem->pair[ict].cont.src[ipart].psirest[icont][BMM_DEM_END_HEAD], M_2PI);

        double zetai = dem->part.r[ipart] * dpsii;
        double zetaj = dem->part.r[jpart] * dpsij;
        double zeta = zetai + zetaj;

        double const proj = bmm_geom2d_dot(xtangij, vdiffij) / d;
        double const dzetai = ri * (dem->part.omega[ipart] - proj);
        double const dzetaj = rj * (dem->part.omega[jpart] - proj);
        double const dzeta = dzetai + dzetaj;

        double const dyn = dem->pair[ict].tang.params.cs.mu * $(bmm_abs, double)(fnorm);
        double const stat = dem->pair[ict].tang.params.cs.k * $(bmm_abs, double)(zeta);

        if (dyn <= stat) {
          // Nothing to do.
        } else {
          e += (1.0 / 2.0) * dem->pair[ict].tang.params.cs.k *
            $(bmm_power, double)(zeta, 2);
        }
      }

      break;
    case BMM_DEM_TANG_BEAM:
      {
        double const lambdaij = bmm_geom2d_dir(xdiffij);
        double const lambdaji = $(bmm_swrap, double)(lambdaij + M_PI, M_2PI);
        double const phii = dem->part.phi[ipart];
        double const phij = dem->part.phi[jpart];
        double const psii = phii - lambdaij;
        double const psij = phij - lambdaji;

        double const dpsii = $(bmm_swrap, double)(psii - dem->pair[ict].cont.src[ipart].psirest[icont][BMM_DEM_END_TAIL], M_2PI);
        double const dpsij = $(bmm_swrap, double)(psij - dem->pair[ict].cont.src[ipart].psirest[icont][BMM_DEM_END_HEAD], M_2PI);

        double zetai = dem->part.r[ipart] * dpsii;
        double zetaj = dem->part.r[jpart] * dpsij;

        e += (1.0 / 2.0) * dem->pair[ict].tang.params.beam.k * $(bmm_power, double)(zetai, 2);
        e += (1.0 / 2.0) * dem->pair[ict].tang.params.beam.k * $(bmm_power, double)(zetaj, 2);
      }

      break;
  }

  return e;
}

// Conservative part only.
__attribute__ ((__nonnull__, __pure__))
double bmm_dem_est_econt(struct bmm_dem const *const dem,
    enum bmm_dem_ct const ict) {
  double e = 0.0;

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    for (size_t icont = 0; icont < dem->pair[ict].cont.src[ipart].n; ++icont) {
      size_t const jpart = dem->pair[ict].cont.src[ipart].itgt[icont];

      e += bmm_dem_est_econt_one(dem, ict, ipart, icont, jpart);
    }

  return e;
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

  double const lambdaij = bmm_geom2d_dir(xdiffij);
  double const lambdaji = $(bmm_swrap, double)(lambdaij + M_PI, M_2PI);
  double const phii = dem->part.phi[ipart];
  double const phij = dem->part.phi[jpart];
  double const psii = phii - lambdaij;
  double const psij = phij - lambdaji;

  dem->pair[ict].cont.src[ipart].strength[icont] = bmm_random_get(dem->rng,
      dem->opts.part.strnew);

  dem->pair[ict].cont.src[ipart].drest[icont] = d;
  dem->pair[ict].cont.src[ipart].itgt[icont] = jpart;

  dem->pair[ict].cont.src[ipart].psirest[icont][BMM_DEM_END_TAIL] = psii;
  dem->pair[ict].cont.src[ipart].psirest[icont][BMM_DEM_END_HEAD] = psij;

  // Not used.
  dem->pair[ict].cont.src[ipart].tfat[icont] = 6;

  double const e = bmm_dem_est_econt_one(dem, ict, ipart, icont, jpart);
  dem->est.ebond += e;
  if (ict == BMM_DEM_CT_WEAK)
    dem->est.ewcont += e;
  else
    dem->est.escont += e;

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
    dem->pair[ict].cont.src[ipart].psirest[icont][iend] = dem->pair[ict].cont.src[ipart].psirest[jcont][iend];

  dem->pair[ict].cont.src[ipart].strength[icont] = dem->pair[ict].cont.src[ipart].strength[jcont];

  dem->pair[ict].cont.src[ipart].tfat[icont] = dem->pair[ict].cont.src[ipart].tfat[jcont];
}

void bmm_dem_remcont_unsafe(struct bmm_dem *const dem,
    enum bmm_dem_ct const ict, size_t const ipart, size_t const icont, size_t const jpart) {
  // fprintf(stderr, "Remove contact %zu out of %zu from %zu to %zu.\n", icont, dem->pair[ict].cont.src[ipart].n, ipart, jpart);

  --dem->pair[ict].cont.src[ipart].n;

  size_t const jcont = dem->pair[ict].cont.src[ipart].n;

  if (jcont != icont)
    bmm_dem_copycont(dem, ict, ipart, icont, jcont);

  // TODO Really?
  // dem->cache.stale = true;
}

void bmm_dem_remcont(struct bmm_dem *const dem,
    enum bmm_dem_ct const ict, size_t apart, size_t bpart, size_t const icont) {
  // TODO Does this ever fail?
  // dynamic_assert(ipart < jpart, "Yes");
  size_t const ipart = $(bmm_min, size_t)(apart, bpart);
  size_t const jpart = $(bmm_max, size_t)(apart, bpart);

  double const e = bmm_dem_est_econt_one(dem, ict, ipart, icont, jpart);
  dem->est.eyieldis += e;
  if (ict == BMM_DEM_CT_WEAK)
    dem->est.ewcont -= e;
  else
    dem->est.escont -= e;

  bmm_dem_remcont_unsafe(dem, ict, ipart, icont, jpart);
}

// I broke this to only work with KV and BEAM.
bool bmm_dem_yield_pair(struct bmm_dem *const dem,
    size_t const ipart, size_t const jpart, size_t const icont) {
  if (dem->part.role[ipart] != BMM_DEM_ROLE_FREE ||
      dem->part.role[jpart] != BMM_DEM_ROLE_FREE)
    return false;

  size_t const ict = BMM_DEM_CT_STRONG;

  double xdiffij[BMM_NDIM];
  bmm_geom2d_cpdiff(xdiffij, dem->part.x[jpart], dem->part.x[ipart],
      dem->opts.box.x, dem->opts.box.per);

  double const ri = dem->part.r[ipart];
  double const rj = dem->part.r[jpart];
  // double const r = ri + rj;
  double const r = dem->pair[ict].cont.src[ipart].drest[icont];
  double const r2 = $(bmm_power, double)(r, 2);

  double const d2 = bmm_geom2d_norm2(xdiffij);
  double const d = sqrt(d2);

  double xnormij[BMM_NDIM];
  bmm_geom2d_scale(xnormij, xdiffij, 1.0 / d);

  double xtangij[BMM_NDIM];
  bmm_geom2d_rperp(xtangij, xnormij);

  double vdiffij[BMM_NDIM];
  bmm_geom2d_diff(vdiffij, dem->part.v[jpart], dem->part.v[ipart]);

  double const xi = r - d;
  double const vnormij = -bmm_geom2d_dot(vdiffij, xnormij);
  double const reff = $(bmm_resum2, double)(ri, rj);
  double const dt = dem->opts.script.dt[dem->script.i];

  double const fnorm = dem->pair[ict].norm.params.dashpot.k * xi +
    dem->pair[ict].norm.params.dashpot.gamma * vnormij;

  double const lambdaij = bmm_geom2d_dir(xdiffij);
  double const lambdaji = $(bmm_swrap, double)(lambdaij + M_PI, M_2PI);
  double const phii = dem->part.phi[ipart];
  double const phij = dem->part.phi[jpart];
  double const psii = phii - lambdaij;
  double const psij = phij - lambdaji;

  double const dpsii = $(bmm_swrap, double)(psii - dem->pair[ict].cont.src[ipart].psirest[icont][BMM_DEM_END_TAIL], M_2PI);
  double const dpsij = $(bmm_swrap, double)(psij - dem->pair[ict].cont.src[ipart].psirest[icont][BMM_DEM_END_HEAD], M_2PI);

  double zetai = dem->part.r[ipart] * dpsii;
  double zetaj = dem->part.r[jpart] * dpsij;

  double const proj = bmm_geom2d_dot(xtangij, vdiffij) / d;
  double const dzetai = ri * (dem->part.omega[ipart] - proj);
  double const dzetaj = rj * (dem->part.omega[jpart] - proj);

  double const ftangconsi = dem->pair[ict].tang.params.beam.k * zetai;
  double const ftangconsj = dem->pair[ict].tang.params.beam.k * zetaj;

  double const ftangdissi = dem->pair[ict].tang.params.beam.dk * dzetai;
  double const ftangdissj = dem->pair[ict].tang.params.beam.dk * dzetaj;

  double const ftangi = ftangconsi + ftangdissi;
  double const ftangj = ftangconsj + ftangdissj;

  double const taui = ri * ftangi;
  double const tauj = rj * ftangj;

  double const ftang = (taui + tauj) / d;

  double const fnormij = fnorm;
  double const ftangij = $(bmm_abs, double)(ftang);

  double const aij = M_PI * $(bmm_power, double)(
      dem->pair[ict].cont.src[ipart].strength[icont] *
      $(bmm_min, double)(dem->part.r[ipart], dem->part.r[jpart]), 2);

  double const sigmanormij = fnormij / aij;
  double const sigmatangij = ftangij / aij;

  double const sigmacrit = fnormij < 0.0 ?
    dem->yield.params.ze.sigmacrit : dem->yield.params.ze.sigmacritt;
  double const taucrit = fnormij < 0.0 ?
    dem->yield.params.ze.taucrit : dem->yield.params.ze.taucritt;

  switch (dem->yield.tag) {
    case BMM_DEM_YIELD_ZE:
      if ($(bmm_power, double)(sigmanormij / sigmacrit, 2) +
          $(bmm_power, double)(sigmatangij / taucrit, 2) > 1.0) {
        bmm_dem_remcont(dem, ict, ipart, jpart, icont);

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

  bool weak = true;
  size_t const icont = bmm_dem_search_cont(dem, BMM_DEM_CT_STRONG, ipart, jpart);
  if (icont != SIZE_MAX)
    weak = bmm_dem_yield_pair(dem, ipart, jpart, icont);

  if (weak) {
    double xdiffij[BMM_NDIM];
    bmm_geom2d_cpdiff(xdiffij, dem->part.x[jpart], dem->part.x[ipart],
        dem->opts.box.x, dem->opts.box.per);

    double const d2 = bmm_geom2d_norm2(xdiffij);
    double const r2 = $(bmm_power, double)(dem->part.r[ipart] + dem->part.r[jpart], 2);
    bool const overlap = d2 < r2;

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
  dem->part.jred[ipart] = bmm_geom_ballprmoi(3);

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

  switch (dem->integ.tag) {
    case BMM_DEM_INTEG_BEEMAN:
      for (size_t idim = 0; idim < BMM_NDIM; ++idim)
        dem->integ.params.beeman.ao[ipart][idim] = 0.0;

      dem->integ.params.beeman.alphao[ipart] = 0.0;

      break;
  }

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

  for (enum bmm_dem_ct ict = 0; ict < BMM_NCT; ++ict) {
    dem->pair[ict].cont.src[ipart].n = dem->pair[ict].cont.src[jpart].n;

    for (size_t icont = 0; icont < dem->pair[ict].cont.src[jpart].n; ++icont)
      dem->pair[ict].cont.src[ipart].itgt[icont] = dem->pair[ict].cont.src[jpart].itgt[icont];

    for (size_t icont = 0; icont < dem->pair[ict].cont.src[jpart].n; ++icont)
      dem->pair[ict].cont.src[ipart].drest[icont] = dem->pair[ict].cont.src[jpart].drest[icont];

    for (size_t icont = 0; icont < dem->pair[ict].cont.src[jpart].n; ++icont)
      for (size_t iend = 0; iend < BMM_NEND; ++iend)
        dem->pair[ict].cont.src[ipart].psirest[icont][iend] = dem->pair[ict].cont.src[jpart].psirest[icont][iend];

    for (size_t icont = 0; icont < dem->pair[ict].cont.src[jpart].n; ++icont)
      dem->pair[ict].cont.src[ipart].strength[icont] = dem->pair[ict].cont.src[jpart].strength[icont];

    for (size_t icont = 0; icont < dem->pair[ict].cont.src[jpart].n; ++icont)
      dem->pair[ict].cont.src[ipart].tfat[icont] = dem->pair[ict].cont.src[jpart].tfat[icont];
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

  double const f = -3.0 * M_2PI * dem->amb.params.creeping.eta *
    dem->part.r[ipart] * v;

  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    dem->part.f[ipart][idim] += f * vunit[idim];

  double const tau = -4.0 * M_2PI * dem->amb.params.creeping.eta *
    $(bmm_power, double)(dem->part.r[ipart], 3);

  dem->part.tau[ipart] += tau * dem->part.omega[ipart];

  double const dt = dem->opts.script.dt[dem->script.i];
  dem->est.eambdis += fabs(f * v * dt);
  dem->est.eambdis += fabs((tau / dem->part.r[ipart]) * v * dt);
}

// TODO Flatten these to allow loop-invariant code motion.
void bmm_dem_force_ambient(struct bmm_dem *const dem, size_t const ipart) {
  switch (dem->amb.tag) {
    case BMM_DEM_AMB_FAXEN:
      bmm_dem_force_creeping(dem, ipart);

      break;
  }
}

void bmm_dem_force_unified(struct bmm_dem *const dem,
    size_t const ipart, size_t const jpart, size_t const icont,
    enum bmm_dem_ct const ict) {
  double xdiffij[BMM_NDIM];
  bmm_geom2d_cpdiff(xdiffij, dem->part.x[jpart], dem->part.x[ipart],
      dem->opts.box.x, dem->opts.box.per);

  double const d2 = bmm_geom2d_norm2(xdiffij);
  if (d2 == 0.0)
    return;

  double const ri = dem->part.r[ipart];
  double const rj = dem->part.r[jpart];
  double const r = ri + rj;
  double const r2 = $(bmm_power, double)(r, 2);
  if (!dem->pair[ict].cohesive)
    if (d2 > r2)
      return;

  double const d = sqrt(d2);

  double xnormij[BMM_NDIM];
  bmm_geom2d_scale(xnormij, xdiffij, 1.0 / d);

  double xtangij[BMM_NDIM];
  bmm_geom2d_rperp(xtangij, xnormij);

  double vdiffij[BMM_NDIM];
  bmm_geom2d_diff(vdiffij, dem->part.v[jpart], dem->part.v[ipart]);

  double dxnorm = 0.0;
  double dxtang = 0.0;

  // Normal forces first.

  double fnorm = 0.0;
  double fnormcons = 0.0;
  double fnormdiss = 0.0;

  {
    double const xi = r - d;
    double const vnormij = -bmm_geom2d_dot(vdiffij, xnormij);
    double const reff = $(bmm_resum2, double)(ri, rj);
    double const dt = dem->opts.script.dt[dem->script.i];
    dxnorm = vnormij * dt;

    switch (dem->pair[ict].norm.tag) {
      case BMM_DEM_NORM_KV:
        if (dem->pair[ict].cohesive) {
          // double const r = ri + rj;
          // Maybe.
          double const r = dem->pair[ict].cont.src[ipart].drest[icont];
          double const xi = r - d;

          fnorm = dem->pair[ict].norm.params.dashpot.k * xi +
            dem->pair[ict].norm.params.dashpot.gamma * vnormij;

          fnormcons = dem->pair[ict].norm.params.dashpot.k * xi;
          fnormdiss = dem->pair[ict].norm.params.dashpot.gamma * vnormij;
        } else {
          fnorm = $(bmm_max, double)(0.0,
              dem->pair[ict].norm.params.dashpot.k * xi +
              dem->pair[ict].norm.params.dashpot.gamma * vnormij);

          if (0.0 > dem->pair[ict].norm.params.dashpot.k * xi +
              dem->pair[ict].norm.params.dashpot.gamma * vnormij) {
            fnormcons = 0.0;
            fnormdiss = 0.0;
          } else {
            fnormcons = dem->pair[ict].norm.params.dashpot.k * xi;
            fnormdiss = dem->pair[ict].norm.params.dashpot.gamma * vnormij;
          }
        }

        break;
      case BMM_DEM_NORM_BSHP:
        {
          double const more = dem->est.bshpp * sqrt(reff * xi);

          if (dem->pair[ict].cohesive) {
            fnorm = more * (xi + dem->pair[ict].norm.params.viscoel.a * vnormij);

            fnormcons = more * xi;
            fnormdiss = more * dem->pair[ict].norm.params.viscoel.a * vnormij;
          } else {
            fnorm = $(bmm_max, double)(0.0,
                more * (xi + dem->pair[ict].norm.params.viscoel.a * vnormij));

            if (0.0 > more * (xi + dem->pair[ict].norm.params.viscoel.a * vnormij)) {
              fnormcons = 0.0;
              fnormdiss = 0.0;
            } else {
              fnormcons = more * xi;
              fnormdiss = more * dem->pair[ict].norm.params.viscoel.a * vnormij;
            }
          }
        }

        break;
    }

    if (ict == BMM_DEM_CT_WEAK) {
      dem->est.ewcont += dxnorm * fnormcons;
      dem->est.ewcontdis += $(bmm_abs, double)(dxnorm * fnormdiss);
    } else {
      dem->est.escont += dxnorm * fnormcons;
      dem->est.escontdis += $(bmm_abs, double)(dxnorm * fnormdiss);
    }
  }

  double fnormij[BMM_NDIM];
  bmm_geom2d_scale(fnormij, xnormij, -fnorm);

  bmm_geom2d_addto(dem->part.f[ipart], fnormij);
  bmm_geom2d_diffto(dem->part.f[jpart], fnormij);

  // Tangential forces second.

  double ftang = 0.0;
  double taui = 0.0;
  double tauj = 0.0;
  double ftangcons = 0.0;
  double ftangdiss = 0.0;

  {
    double const dt = dem->opts.script.dt[dem->script.i];

    double const lambdaij = bmm_geom2d_dir(xdiffij);
    double const lambdaji = $(bmm_swrap, double)(lambdaij + M_PI, M_2PI);
    double const phii = dem->part.phi[ipart];
    double const phij = dem->part.phi[jpart];
    double const psii = phii - lambdaij;
    double const psij = phij - lambdaji;

    double const dpsii = $(bmm_swrap, double)(psii - dem->pair[ict].cont.src[ipart].psirest[icont][BMM_DEM_END_TAIL], M_2PI);
    double const dpsij = $(bmm_swrap, double)(psij - dem->pair[ict].cont.src[ipart].psirest[icont][BMM_DEM_END_HEAD], M_2PI);

    double zetai = dem->part.r[ipart] * dpsii;
    double zetaj = dem->part.r[jpart] * dpsij;
    double zeta = zetai + zetaj;

    double const proj = bmm_geom2d_dot(xtangij, vdiffij) / d;
    double const dzetai = ri * (dem->part.omega[ipart] - proj);
    double const dzetaj = rj * (dem->part.omega[jpart] - proj);
    double const dzeta = dzetai + dzetaj;

    dxtang = dzeta * dt;

    switch (dem->pair[ict].tang.tag) {
      case BMM_DEM_TANG_HW:
        {
          double const dyn = dem->pair[ict].tang.params.hw.mu * $(bmm_abs, double)(fnorm);
          double const shear = dem->pair[ict].tang.params.hw.gamma * $(bmm_abs, double)(dzeta);

          ftang = copysign($(bmm_min, double)(dyn, shear), dzeta);
          if (dyn <= shear)
            ++dem->est.hwmu;
          else
            ++dem->est.hwgamma;

          ftangcons = 0.0;
          ftangdiss = ftang;

          taui = ri * ftang;
          tauj = rj * ftang;
        }
        break;
      case BMM_DEM_TANG_CS:
        {
          double const dyn = dem->pair[ict].tang.params.cs.mu * $(bmm_abs, double)(fnorm);
          double const stat = dem->pair[ict].tang.params.cs.k * $(bmm_abs, double)(zeta);

          ftang = dyn <= stat ? copysign(dyn, dzeta) : copysign(stat, zeta);

          if (dyn <= stat) {
            ftangcons = 0.0;
            ftangdiss = copysign(dyn, dzeta);

            ++dem->est.csmu;
          } else {
            ftangcons = copysign(stat, zeta);
            ftangdiss = 0.0;

            ++dem->est.csk;
          }

          taui = ri * ftang;
          tauj = rj * ftang;
        }

        break;
      case BMM_DEM_TANG_BEAM:
        {
          double const d = bmm_geom2d_norm(xdiffij);
          double const r = ri + rj;
          double const r2 = $(bmm_power, double)(r, 2);

          double const lambdaij = bmm_geom2d_dir(xdiffij);
          double const lambdaji = $(bmm_swrap, double)(lambdaij + M_PI, M_2PI);
          double const phii = dem->part.phi[ipart];
          double const phij = dem->part.phi[jpart];
          double const psii = phii - lambdaij;
          double const psij = phij - lambdaji;

          double const dpsii = $(bmm_swrap, double)(psii - dem->pair[ict].cont.src[ipart].psirest[icont][BMM_DEM_END_TAIL], M_2PI);
          double const dpsij = $(bmm_swrap, double)(psij - dem->pair[ict].cont.src[ipart].psirest[icont][BMM_DEM_END_HEAD], M_2PI);

          double zetai = dem->part.r[ipart] * dpsii;
          double zetaj = dem->part.r[jpart] * dpsij;

          double const proj = bmm_geom2d_dot(xtangij, vdiffij) / d;
          double const dzetai = ri * (dem->part.omega[ipart] - proj);
          double const dzetaj = rj * (dem->part.omega[jpart] - proj);

          double const ftangconsi = dem->pair[ict].tang.params.beam.k * zetai;
          double const ftangconsj = dem->pair[ict].tang.params.beam.k * zetaj;

          double const ftangdissi = dem->pair[ict].tang.params.beam.dk * dzetai;
          double const ftangdissj = dem->pair[ict].tang.params.beam.dk * dzetaj;

          double const ftangi = ftangconsi + ftangdissi;
          double const ftangj = ftangconsj + ftangdissj;

          taui = ri * ftangi;
          tauj = rj * ftangj;

          ftang = (taui + tauj) / d;

          if (ict == BMM_DEM_CT_WEAK) {
            dem->est.ewcont += (dzetai * dt) * ftangconsi;
            dem->est.ewcont += (dzetaj * dt) * ftangconsj;
            dem->est.ewcontdis += $(bmm_abs, double)((dzetai * dt) * ftangdissi);
            dem->est.ewcontdis += $(bmm_abs, double)((dzetaj * dt) * ftangdissj);
          } else {
            dem->est.escont += (dzetai * dt) * ftangconsi;
            dem->est.escont += (dzetaj * dt) * ftangconsj;
            dem->est.escontdis += $(bmm_abs, double)((dzetai * dt) * ftangdissi);
            dem->est.escontdis += $(bmm_abs, double)((dzetaj * dt) * ftangdissj);
          }
        }

        break;
    }

    if (ict == BMM_DEM_CT_WEAK) {
      dem->est.ewcont += dxtang * ftangcons;
      dem->est.ewcontdis += $(bmm_abs, double)(dxtang * ftangdiss);
    } else {
      dem->est.escont += dxtang * ftangcons;
      dem->est.escontdis += $(bmm_abs, double)(dxtang * ftangdiss);
    }
  }

  double ftangij[BMM_NDIM];
  bmm_geom2d_scale(ftangij, xtangij, -ftang);

  bmm_geom2d_addto(dem->part.f[ipart], ftangij);
  bmm_geom2d_diffto(dem->part.f[jpart], ftangij);

  dem->part.tau[ipart] -= taui;
  dem->part.tau[jpart] -= tauj;
}

void bmm_dem_force_external(struct bmm_dem *const dem, size_t const ipart) {
  double const dt = dem->opts.script.dt[dem->script.i];

  switch (dem->ext.tag) {
    case BMM_DEM_EXT_HARM:
      dem->part.f[ipart][1] += dem->ext.params.harm.k *
        (dem->opts.box.x[1] / 2.0 - dem->part.x[ipart][1]);

      break;
    case BMM_DEM_EXT_GRAVY:
      {
        dem->part.f[ipart][1] += dem->part.m[ipart] * dem->ext.params.gravy.g;

        double const dy = dem->part.v[ipart][1] * dt;
      }

      break;
    case BMM_DEM_EXT_DRIVE:
      switch (dem->part.role[ipart]) {
        case BMM_DEM_ROLE_DRIVEN:
          {
            double f[BMM_NDIM];
            for (size_t idim = 0; idim < BMM_NDIM; ++idim)
              f[idim] = dem->script.state.crunch.fdrive[idim] /
                (double) dem->script.state.crunch.ndrive;

            double const dx = dem->part.v[ipart][0] * dt;
            double const dy = dem->part.v[ipart][1] * dt;
            dem->est.edrivtang += dx * f[0];
            dem->est.edrivnorm += dy * f[1];

            for (size_t idim = 0; idim < BMM_NDIM; ++idim)
              dem->part.f[ipart][idim] += f[idim];
          }

          break;
      }

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

  // This goes for the previous frame, so this is not exactly the right spot.
  // dem->est.eambdis = 0.0;
  dem->est.epotext_d = 0.0;
  dem->est.eklin_d = 0.0;
  dem->est.ekrot_d = 0.0;
  dem->est.ewcont_d = 0.0;
  dem->est.escont_d = 0.0;
  // dem->est.ewcont = 0.0;
  // dem->est.escont = 0.0;
  // dem->est.edrivtang = 0.0;
  // dem->est.edrivnorm = 0.0;
  // dem->est.ebond = 0.0;
  // dem->est.eyieldis = 0.0;
  // dem->est.ewcontdis = 0.0;
  // dem->est.escontdis = 0.0;
  dem->est.fback[0] = 0.0;
  dem->est.fback[1] = 0.0;
  dem->est.mueff = (double) NAN;
  dem->est.mueffb = (double) NAN;
  dem->est.vdriv[0] = (double) NAN;
  dem->est.vdriv[1] = (double) NAN;
  dem->est.bshpp = (2.0 / 3.0) * (dem->opts.part.ycomp /
      (1.0 - $(bmm_power, double)(dem->opts.part.nu, 2)));


  dem->est.epotext_d += bmm_dem_est_epotext(dem);
  dem->est.eklin_d += bmm_dem_est_eklin(dem);
  dem->est.ekrot_d += bmm_dem_est_ekrot(dem);
  dem->est.ewcont_d += bmm_dem_est_econt(dem, BMM_DEM_CT_WEAK);
  dem->est.escont_d += bmm_dem_est_econt(dem, BMM_DEM_CT_STRONG);

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    bmm_dem_force_ambient(dem, ipart);

  for (enum bmm_dem_ct ict = 0; ict < BMM_NCT; ++ict)
    for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
      for (size_t icont = 0; icont < dem->pair[ict].cont.src[ipart].n; ++icont) {
        size_t const jpart = dem->pair[ict].cont.src[ipart].itgt[icont];

        bmm_dem_force_unified(dem, ipart, jpart, icont, ict);
      }

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    if (dem->part.role[ipart] == BMM_DEM_ROLE_DRIVEN)
      for (size_t idim = 0; idim < BMM_NDIM; ++idim)
        dem->est.fback[idim] += dem->part.f[ipart][idim];

  if (dem->script.state.crunch.fdrive[1] != 0.0)
    dem->est.mueff = $(bmm_abs, double)(dem->script.state.crunch.fdrive[0] /
        dem->script.state.crunch.fdrive[1]);

  if (dem->est.fback[1] != 0.0)
    dem->est.mueffb = $(bmm_abs, double)(dem->est.fback[0] /
        dem->est.fback[1]);

  bmm_dem_est_vdc(dem->est.vdriv, dem);

  dem->script.state.crunch.fdrive[0] +=
    copysign(dem->opts.script.params[dem->script.i].crunch.fadjust[0] * dt,
        dem->opts.script.params[dem->script.i].crunch.v - dem->est.vdriv[0]);

  dem->script.state.crunch.fdrive[1] +=
    copysign(dem->opts.script.params[dem->script.i].crunch.fadjust[1] * dt,
        dem->est.fback[1] -
        dem->opts.script.params[dem->script.i].crunch.p *
        dem->opts.box.x[0] *
        2.0 * bmm_ival_midpoint(dem->opts.part.rnew));

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    bmm_dem_force_external(dem, ipart);
}

void bmm_dem_accel(struct bmm_dem *const dem) {
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    if (dem->part.role[ipart] != BMM_DEM_ROLE_FIXED) {
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

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    for (size_t idim = 0; idim < BMM_NDIM; ++idim) {
      dem->part.x[ipart][idim] += dem->part.v[ipart][idim] * dt;
      dem->part.v[ipart][idim] += dem->part.a[ipart][idim] * dt;

      if (dem->opts.box.per[idim])
        dem->part.x[ipart][idim] = $(bmm_uwrap, double)(dem->part.x[ipart][idim],
            dem->opts.box.x[idim]);
    }

    dem->part.phi[ipart] += dem->part.omega[ipart] * dt;
    dem->part.omega[ipart] += dem->part.alpha[ipart] * dt;
  }
}

void bmm_dem_integ_taylor(struct bmm_dem *const dem) {
  double const dt = dem->opts.script.dt[dem->script.i];

  if (dt == 0.0)
    return;

  double const dt2 = $(bmm_power, double)(dt, 2);

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    for (size_t idim = 0; idim < BMM_NDIM; ++idim) {
      dem->part.x[ipart][idim] += dem->part.v[ipart][idim] * dt +
        (1.0 / 2.0) * dem->part.a[ipart][idim] * dt2;
      dem->part.v[ipart][idim] += dem->part.a[ipart][idim] * dt;

      if (dem->opts.box.per[idim])
        dem->part.x[ipart][idim] = $(bmm_uwrap, double)(dem->part.x[ipart][idim],
            dem->opts.box.x[idim]);
    }

    dem->part.phi[ipart] += dem->part.omega[ipart] * dt +
      (1.0 / 2.0) * dem->part.alpha[ipart] * dt2;
    dem->part.omega[ipart] += dem->part.alpha[ipart] * dt;
  }
}

void bmm_dem_integ_vel(struct bmm_dem *const dem) {
  double const dt = dem->opts.script.dt[dem->script.i];

  if (dt == 0.0)
    return;

  double const dt2 = $(bmm_power, double)(dt, 2);

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    for (size_t idim = 0; idim < BMM_NDIM; ++idim) {
      dem->integ.params.velvet.ao[ipart][idim] = dem->part.a[ipart][idim];

      dem->part.x[ipart][idim] += dem->part.v[ipart][idim] * dt +
        (1.0 / 2.0) * dem->part.a[ipart][idim] * dt2;

      if (dem->opts.box.per[idim])
        dem->part.x[ipart][idim] = $(bmm_uwrap, double)(dem->part.x[ipart][idim],
            dem->opts.box.x[idim]);
    }

    dem->integ.params.velvet.alphao[ipart] = dem->part.alpha[ipart];

    dem->part.phi[ipart] += dem->part.omega[ipart] * dt +
      (1.0 / 2.0) * dem->part.alpha[ipart] * dt2;
  }
}

void bmm_dem_integ_vet(struct bmm_dem *const dem) {
  double const dt = dem->opts.script.dt[dem->script.i];

  if (dt == 0.0)
    return;

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    for (size_t idim = 0; idim < BMM_NDIM; ++idim)
      dem->part.v[ipart][idim] += (1.0 / 2.0) *
        (dem->part.a[ipart][idim] + dem->integ.params.velvet.ao[ipart][idim]) * dt;

    dem->part.omega[ipart] += (1.0 / 2.0) *
      (dem->part.alpha[ipart] + dem->integ.params.velvet.alphao[ipart]) * dt;
  }
}

void bmm_dem_integ_bee(struct bmm_dem *const dem) {
  double const dt = dem->opts.script.dt[dem->script.i];

  if (dt == 0.0)
    return;

  double const dt2 = $(bmm_power, double)(dt, 2);

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    for (size_t idim = 0; idim < BMM_NDIM; ++idim) {
      dem->integ.params.beeman.xo[ipart][idim] = dem->part.x[ipart][idim];
      dem->integ.params.beeman.vo[ipart][idim] = dem->part.v[ipart][idim];
      dem->integ.params.beeman.aoo[ipart][idim] = dem->integ.params.beeman.ao[ipart][idim];
      dem->integ.params.beeman.ao[ipart][idim] = dem->part.a[ipart][idim];

      dem->part.x[ipart][idim] = dem->part.x[ipart][idim] +
        dem->part.v[ipart][idim] * dt +
        (1.0 / 2.0) * (1.0 / 3.0) * (4.0 * dem->part.a[ipart][idim] -
            dem->integ.params.beeman.aoo[ipart][idim]) * dt2;

      if (dem->opts.box.per[idim])
        dem->part.x[ipart][idim] = $(bmm_uwrap, double)(dem->part.x[ipart][idim],
            dem->opts.box.x[idim]);

      dem->part.v[ipart][idim] = dem->part.v[ipart][idim] +
        (1.0 / 2.0) * (3.0 * dem->part.a[ipart][idim] -
            dem->integ.params.beeman.aoo[ipart][idim]) * dt;
    }

    dem->integ.params.beeman.phio[ipart] = dem->part.phi[ipart];
    dem->integ.params.beeman.omegao[ipart] = dem->part.omega[ipart];
    dem->integ.params.beeman.alphaoo[ipart] = dem->integ.params.beeman.alphao[ipart];
    dem->integ.params.beeman.alphao[ipart] = dem->part.alpha[ipart];

    dem->part.phi[ipart] = dem->part.phi[ipart] +
      dem->part.omega[ipart] * dt +
      (1.0 / 2.0) * (1.0 / 3.0) * (4.0 * dem->part.alpha[ipart] -
          dem->integ.params.beeman.alphaoo[ipart]) * dt2;

    dem->part.omega[ipart] = dem->part.omega[ipart] +
      (1.0 / 2.0) * (3.0 * dem->part.alpha[ipart] -
          dem->integ.params.beeman.alphaoo[ipart]) * dt;
  }
}

void bmm_dem_integ_man(struct bmm_dem *const dem) {
  double const dt = dem->opts.script.dt[dem->script.i];

  if (dt == 0.0)
    return;

  double const dt2 = $(bmm_power, double)(dt, 2);

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    for (size_t idim = 0; idim < BMM_NDIM; ++idim) {
      dem->part.x[ipart][idim] = dem->integ.params.beeman.xo[ipart][idim] +
        dem->integ.params.beeman.vo[ipart][idim] * dt +
        (1.0 / 2.0) * (1.0 / 3.0) * (dem->part.a[ipart][idim] +
            2.0 * dem->integ.params.beeman.ao[ipart][idim]) * dt2;

      if (dem->opts.box.per[idim])
        dem->part.x[ipart][idim] = $(bmm_uwrap, double)(dem->part.x[ipart][idim],
            dem->opts.box.x[idim]);

      dem->part.v[ipart][idim] = dem->integ.params.beeman.vo[ipart][idim] +
        (1.0 / 6.0) * (2.0 * dem->part.a[ipart][idim] +
            5.0 * dem->integ.params.beeman.ao[ipart][idim] -
            dem->integ.params.beeman.aoo[ipart][idim]) * dt;
    }

    dem->part.phi[ipart] = dem->integ.params.beeman.phio[ipart] +
      dem->integ.params.beeman.omegao[ipart] * dt +
      (1.0 / 2.0) * (1.0 / 3.0) * (dem->part.alpha[ipart] +
          2.0 * dem->integ.params.beeman.alphao[ipart]) * dt2;

    dem->part.omega[ipart] = dem->integ.params.beeman.omegao[ipart] +
      (1.0 / 6.0) * (2.0 * dem->part.alpha[ipart] +
          5.0 * dem->integ.params.beeman.alphao[ipart] -
          dem->integ.params.beeman.alphaoo[ipart]) * dt;
  }
}

void bmm_dem_integ_kura(struct bmm_dem *const dem) {
  double const dt = dem->opts.script.dt[dem->script.i];

  if (dt == 0.0)
    return;

  double const dt2 = $(bmm_power, double)(dt, 2);

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    for (size_t idim = 0; idim < BMM_NDIM; ++idim) {
      dem->integ.params.beeman.xo[ipart][idim] = dem->part.x[ipart][idim];
      dem->integ.params.beeman.vo[ipart][idim] = dem->part.v[ipart][idim];
      dem->integ.params.beeman.aoo[ipart][idim] = dem->integ.params.beeman.ao[ipart][idim];
      dem->integ.params.beeman.ao[ipart][idim] = dem->part.a[ipart][idim];

      dem->part.x[ipart][idim] = dem->part.x[ipart][idim] +
        dem->part.v[ipart][idim] * dt +
        (1.0 / 2.0) * (1.0 / 4.0) * (5.0 * dem->part.a[ipart][idim] -
            dem->integ.params.beeman.aoo[ipart][idim]) * dt2;

      if (dem->opts.box.per[idim])
        dem->part.x[ipart][idim] = $(bmm_uwrap, double)(dem->part.x[ipart][idim],
            dem->opts.box.x[idim]);

      dem->part.v[ipart][idim] = dem->part.v[ipart][idim] +
        (1.0 / 2.0) * (3.0 * dem->part.a[ipart][idim] -
            dem->integ.params.beeman.aoo[ipart][idim]) * dt;
    }

    dem->integ.params.beeman.phio[ipart] = dem->part.phi[ipart];
    dem->integ.params.beeman.omegao[ipart] = dem->part.omega[ipart];
    dem->integ.params.beeman.alphaoo[ipart] = dem->integ.params.beeman.alphao[ipart];
    dem->integ.params.beeman.alphao[ipart] = dem->part.alpha[ipart];

    dem->part.phi[ipart] = dem->part.phi[ipart] +
      dem->part.omega[ipart] * dt +
      (1.0 / 2.0) * (1.0 / 3.0) * (4.0 * dem->part.alpha[ipart] -
          dem->integ.params.beeman.alphaoo[ipart]) * dt2;

    dem->part.omega[ipart] = dem->part.omega[ipart] +
      (1.0 / 2.0) * (3.0 * dem->part.alpha[ipart] -
          dem->integ.params.beeman.alphaoo[ipart]) * dt;
  }
}

void bmm_dem_integ_ev(struct bmm_dem *const dem) {
  double const dt = dem->opts.script.dt[dem->script.i];

  if (dt == 0.0)
    return;

  double const dt2 = $(bmm_power, double)(dt, 2);

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    for (size_t idim = 0; idim < BMM_NDIM; ++idim) {
      dem->part.x[ipart][idim] = dem->integ.params.beeman.xo[ipart][idim] +
        dem->integ.params.beeman.vo[ipart][idim] * dt +
        (1.0 / 2.0) * (1.0 / 4.0) * (dem->part.a[ipart][idim] +
            3.0 * dem->integ.params.beeman.ao[ipart][idim]) * dt2;

      if (dem->opts.box.per[idim])
        dem->part.x[ipart][idim] = $(bmm_uwrap, double)(dem->part.x[ipart][idim],
            dem->opts.box.x[idim]);

      dem->part.v[ipart][idim] = dem->integ.params.beeman.vo[ipart][idim] +
        (1.0 / 8.0) * (3.0 * dem->part.a[ipart][idim] +
            6.0 * dem->integ.params.beeman.ao[ipart][idim] -
            dem->integ.params.beeman.aoo[ipart][idim]) * dt;
    }

    dem->part.phi[ipart] = dem->integ.params.beeman.phio[ipart] +
      dem->integ.params.beeman.omegao[ipart] * dt +
      (1.0 / 2.0) * (1.0 / 3.0) * (dem->part.alpha[ipart] +
          2.0 * dem->integ.params.beeman.alphao[ipart]) * dt2;

    dem->part.omega[ipart] = dem->integ.params.beeman.omegao[ipart] +
      (1.0 / 6.0) * (2.0 * dem->part.alpha[ipart] +
          5.0 * dem->integ.params.beeman.alphao[ipart] -
          dem->integ.params.beeman.alphaoo[ipart]) * dt;
  }
}

void bmm_dem_stab(struct bmm_dem *const dem) {
  // We could reset windings here,
  // but that would interfere with beam models,
  // so we relax instead.
}

void bmm_dem_predict(struct bmm_dem *const dem) {
  switch (dem->integ.tag) {
    case BMM_DEM_INTEG_VELVET:
      bmm_dem_integ_vel(dem);

      break;
    case BMM_DEM_INTEG_BEEMAN:
      bmm_dem_integ_bee(dem);

      break;
    case BMM_DEM_INTEG_KURAEV:
      bmm_dem_integ_kura(dem);

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
    case BMM_DEM_INTEG_BEEMAN:
      bmm_dem_integ_man(dem);

      break;
    case BMM_DEM_INTEG_KURAEV:
      bmm_dem_integ_ev(dem);

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

  if (d2 > r2 * dem->bond.ccrcont)
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
    for (enum bmm_dem_ct ict = 0; ict < BMM_NCT; ++ict)
      dem->pair[ict].cont.src[ipart].n = 0;
}

__attribute__ ((__nonnull__))
double bmm_dem_bulkmod(struct bmm_dem const *const dem) {
  return dem->opts.part.ycomp / (3.0 * (1.0 - 2.0 * dem->opts.part.nu));
}

__attribute__ ((__nonnull__))
double bmm_dem_shearmod(struct bmm_dem const *const dem) {
  return dem->opts.part.ycomp / (2.0 * (1.0 + dem->opts.part.nu));
}

__attribute__ ((__nonnull__))
double bmm_dem_pwavec(struct bmm_dem const *const dem) {
  return sqrt((bmm_dem_bulkmod(dem) +
        (4.0 / 3.0) * bmm_dem_shearmod(dem)) / dem->opts.part.rho);
}

__attribute__ ((__nonnull__))
double bmm_dem_swavec(struct bmm_dem const *const dem) {
  return sqrt(bmm_dem_shearmod(dem) / dem->opts.part.rho);
}

static bool export_s(struct bmm_dem const *const dem) {
  char buf[BUFSIZ];
  (void) snprintf(buf, sizeof buf, "./%s-s.data",
      dem->opts.script.params[dem->script.i].expr.str);

  FILE *const stream = fopen(buf, "w");
  if (stream == NULL) {
    BMM_TLE_STDS();

    return false;
  }

  if (fprintf(stream, "%g\n", dem->est.sk) < 0) {
    BMM_TLE_STDS();

    return false;
  }

  if (fclose(stream) != 0) {
    BMM_TLE_STDS();

    return false;
  }

  return true;
}

static bool export_chi(struct bmm_dem const *const dem) {
  char buf[BUFSIZ];
  (void) snprintf(buf, sizeof buf, "./%s-chi.data",
      dem->opts.script.params[dem->script.i].expr.str);

  FILE *const stream = fopen(buf, "w");
  if (stream == NULL) {
    BMM_TLE_STDS();

    return false;
  }

  if (fprintf(stream, "%g\n", dem->est.chi) < 0) {
    BMM_TLE_STDS();

    return false;
  }

  if (fclose(stream) != 0) {
    BMM_TLE_STDS();

    return false;
  }

  return true;
}

static bool export_x(struct bmm_dem const *const dem) {
  char buf[BUFSIZ];
  (void) snprintf(buf, sizeof buf, "./%s-x.data",
      dem->opts.script.params[dem->script.i].expr.str);

  FILE *const stream = fopen(buf, "w");
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

static bool export_phi(struct bmm_dem const *const dem) {
  char buf[BUFSIZ];
  (void) snprintf(buf, sizeof buf, "./%s-phi.data",
      dem->opts.script.params[dem->script.i].expr.str);

  FILE *const stream = fopen(buf, "w");
  if (stream == NULL) {
    BMM_TLE_STDS();

    return false;
  }

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    if (fprintf(stream, "%g\n", dem->part.phi[ipart]) < 0) {
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
  char buf[BUFSIZ];
  (void) snprintf(buf, sizeof buf, "./%s-r.data",
      dem->opts.script.params[dem->script.i].expr.str);

  FILE *const stream = fopen(buf, "w");
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
  char buf[BUFSIZ];
  (void) snprintf(buf, sizeof buf, "./%s-c.data",
      dem->opts.script.params[dem->script.i].expr.str);

  FILE *const stream = fopen(buf, "w");
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
  char buf[BUFSIZ];
  (void) snprintf(buf, sizeof buf, "./%s-f.data",
      dem->opts.script.params[dem->script.i].expr.str);

  FILE *const stream = fopen(buf, "w");
  if (stream == NULL) {
    BMM_TLE_STDS();

    return false;
  }

  // Bidirectionalization of the strong contact graph.
  // Following great programming conventions, I copied this from below.
  struct agraph *const agraph = malloc(sizeof *agraph);
  if (agraph == NULL) {
    BMM_TLE_STDS();

    return false;
  }
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
  char buf[BUFSIZ];
  (void) snprintf(buf, sizeof buf, "./%s-p.data",
      dem->opts.script.params[dem->script.i].expr.str);

  FILE *const stream = fopen(buf, "w");
  if (stream == NULL) {
    BMM_TLE_STDS();

    return false;
  }

  // Bidirectionalization of the strong contact graph.
  struct agraph *const agraph = malloc(sizeof *agraph);
  if (agraph == NULL) {
    BMM_TLE_STDS();

    return false;
  }
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
  if (faces == NULL) {
    BMM_TLE_STDS();

    return false;
  }
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

    double backdir = (double) NAN;
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

// TODO These procedures do not belong here.

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

  // double g[BMM_TRI(BMM_MPART) * 2];
  size_t nmemb = $(bmm_tri, size_t)(dem->part.n) * 2;
  double *const w = malloc(nmemb * sizeof *w);
  if (w == NULL)
    return false;
  double *const r = malloc(nmemb * sizeof *r);
  if (r == NULL) {
    free(w);
    return false;
  }

// If you desperately want wrong results...
// #define NO_WEIGHTS

  nmemb = 0;
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    if (bmm_dem_inside(dem, ipart))
      for (size_t jpart = 0; jpart < dem->part.n; ++jpart)
        if (jpart != ipart && bmm_dem_inside(dem, jpart)) {
          double const d = bmm_geom2d_cpdist(dem->part.x[ipart],
              dem->part.x[jpart], dem->opts.box.x, dem->opts.box.per);

          double const v0 = bmm_geom_ballsurf(d, 2);
          double const v = bmm_geom2d_shellvol(dem->part.x[ipart], d,
              dem->opts.box.x, dem->opts.box.per);

          r[nmemb] = d;
#ifdef NO_WEIGHTS
          w[nmemb] = 1.0;
#else
          w[nmemb] = d == 0.0 ? 0.0 : v0 / v;
#endif
          ++nmemb;
        }

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

  // We must renormalize to account for cut off kernels.
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

static bool dump_raddist_etc(struct bmm_dem *const dem) {
  FILE *const stream = fopen("raddist.data", "w");
  if (stream == NULL) {
    BMM_TLE_STDS();

    return false;
  }

  size_t nbin = 512;
  double *const r = malloc(nbin * sizeof *r);
  double *const g = malloc(nbin * sizeof *g);
  dynamic_assert(r != NULL && g != NULL, "Allocated");

  double const rmax = dem->opts.box.x[0] / 4.0;

  if (!bmm_dem_est_raddist(r, g, nbin, rmax, dem))
    return false;

  for (size_t ibin = 0; ibin < nbin; ++ibin)
    if (fprintf(stream, "%g %g\n", r[ibin], g[ibin]) < 0) {
      BMM_TLE_STDS();

      return false;
    }

  // Snip.

  double const v = $(bmm_prod, double)(dem->opts.box.x, BMM_NDIM);
  double const rho = (double) dem->part.n / v;
  double const dr = rmax / (double) nbin;
  double a = 0.0;

  for (size_t ibin = 0; ibin < nbin; ++ibin)
    a += (g[ibin] * log(g[ibin] + 1.0e-9) - (g[ibin] - 1.0)) * dr;

  dem->est.sk = -(rho / 2.0) * a;

  // Snap.

  free(g);
  free(r);

  if (fclose(stream) != 0) {
    BMM_TLE_STDS();

    return false;
  }

  return true;
}

void bmm_dem_opts_def(struct bmm_dem_opts *const opts) {
  // This is here just to help Valgrind and cover up my mistakes.
  (void) memset(opts, 0, sizeof *opts);

  opts->verbose = false;

  opts->gross.fric = BMM_DEM_TANG_CS;
  opts->gross.pfac = 0.5;
  opts->gross.hfac = 1.0;
  opts->gross.ds = 0.1;

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

  opts->part.ytens = 1.0;
  opts->part.ycomp = 1.0;
  opts->part.nu = 0.5;
  opts->part.rnew[0] = 1.0;
  opts->part.rnew[1] = 1.0;
  opts->part.strnew[0] = 1.0;
  opts->part.strnew[1] = 1.0;

  opts->script.n = 0;

  opts->comm.dt = 1.0;
  opts->comm.flip = true;
  opts->comm.flop = true;
  opts->comm.flap = true;
  opts->comm.flup = true;

  opts->cache.dcutoff = 1.0 / 5.0;

  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    opts->cache.ncell[idim] = 5;
}

void bmm_dem_def(struct bmm_dem *const dem,
    struct bmm_dem_opts const *const opts) {
  // This is here just to help Valgrind and cover up my mistakes.
  (void) memset(dem, 0, sizeof *dem);

  dem->opts = *opts;

  dem->trap.remask = 0;

  dem->bond.ccrcont = 1.5;

  dem->integ.tag = BMM_DEM_INTEG_TAYLOR;
  dem->integ.tag = BMM_DEM_INTEG_VELVET;
  dem->integ.tag = BMM_DEM_INTEG_BEEMAN;
  dem->integ.tag = BMM_DEM_INTEG_KURAEV;
  dem->integ.tag = BMM_DEM_INTEG_EULER;
  dem->cache.tag = BMM_DEM_CACHE_NONE;
  dem->cache.tag = BMM_DEM_CACHE_NEIGH;
  dem->ext.tag = BMM_DEM_EXT_NONE;
  dem->amb.tag = BMM_DEM_AMB_NONE;
  dem->pair[BMM_DEM_CT_WEAK].cohesive = false;
  dem->pair[BMM_DEM_CT_WEAK].norm.tag = BMM_DEM_NORM_BSHP;
  dem->pair[BMM_DEM_CT_WEAK].tang.tag = BMM_DEM_TANG_CS;
  dem->pair[BMM_DEM_CT_STRONG].cohesive = true;
  dem->pair[BMM_DEM_CT_STRONG].norm.tag = BMM_DEM_NORM_KV;
  dem->pair[BMM_DEM_CT_STRONG].tang.tag = BMM_DEM_TANG_BEAM;
  dem->yield.tag = BMM_DEM_YIELD_ZE;

  switch (dem->amb.tag) {
    case BMM_DEM_AMB_FAXEN:
      dem->amb.params.creeping.eta = 1.0e+3;

      break;
  }

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
      dem->pair[BMM_DEM_CT_WEAK].tang.params.hw.mu = 0.85;

      break;
    case BMM_DEM_TANG_CS:
      dem->pair[BMM_DEM_CT_WEAK].tang.params.cs.k = 1.0e+6;
      dem->pair[BMM_DEM_CT_WEAK].tang.params.cs.mu = 0.85;

      break;
  }

  switch (dem->pair[BMM_DEM_CT_STRONG].norm.tag) {
    case BMM_DEM_NORM_KV:
      dem->pair[BMM_DEM_CT_STRONG].norm.params.dashpot.k = 2.0e+8;
      dem->pair[BMM_DEM_CT_STRONG].norm.params.dashpot.gamma = 4.0e+3;

      break;
    case BMM_DEM_NORM_BSHP:
      dem->pair[BMM_DEM_CT_STRONG].norm.params.viscoel.a = 4.0e-2;

      break;
  }

  switch (dem->pair[BMM_DEM_CT_STRONG].tang.tag) {
    case BMM_DEM_TANG_BEAM:
      dem->pair[BMM_DEM_CT_STRONG].tang.params.beam.k = 8.0e+6;
      dem->pair[BMM_DEM_CT_STRONG].tang.params.beam.dk = 4.0e+2;
      // This is ideal for `triplet` tests.
      // dem->pair[BMM_DEM_CT_STRONG].tang.params.beam.dk = 4.0e+1;
      // dem->pair[BMM_DEM_CT_STRONG].tang.params.beam.dk = 0.0e+1;

      break;
  }

  switch (dem->yield.tag) {
    case BMM_DEM_YIELD_RANKINE:
      dem->yield.params.rankine.sigmacrit = 1.0e+8;

      break;
    case BMM_DEM_YIELD_TRESCA:
      dem->yield.params.tresca.taucrit = 1.0e+8;

      break;
    case BMM_DEM_YIELD_ZE:
      dem->yield.params.ze.sigmacrit = 5.0e+8;
      dem->yield.params.ze.taucrit = 5.0e+8;
      // This is ideal for `beam` tests.
      // dem->yield.params.ze.sigmacrit = 5.0e+7;
      // dem->yield.params.ze.taucrit = 5.0e+7;

      break;
  }

  // Estimators appeared.
  dem->est.eambdis = 0.0;
  dem->est.epotext_d = 0.0;
  dem->est.eklin_d = 0.0;
  dem->est.ekrot_d = 0.0;
  dem->est.ewcont_d = 0.0;
  dem->est.escont_d = 0.0;
  dem->est.ewcont = 0.0;
  dem->est.escont = 0.0;
  dem->est.edrivtang = 0.0;
  dem->est.edrivnorm = 0.0;
  dem->est.ebond = 0.0;
  dem->est.eyieldis = 0.0;
  dem->est.ewcontdis = 0.0;
  dem->est.escontdis = 0.0;
  dem->est.fback[0] = 0.0;
  dem->est.fback[1] = 0.0;
  dem->est.chi = (double) NAN;
  dem->est.mueff = (double) NAN;
  dem->est.mueffb = (double) NAN;
  dem->est.vdriv[0] = (double) NAN;
  dem->est.vdriv[1] = (double) NAN;
  dem->est.hwgamma = 0;
  dem->est.hwmu = 0;
  dem->est.csk = 0;
  dem->est.csmu = 0;

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
static bool bmm_dem_script_set_density(struct bmm_dem *const dem) {
  double const vhc = $(bmm_prod, double)(dem->opts.box.x, BMM_NDIM);

  double v = 0.0;
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    v += bmm_geom_ballvol(dem->part.r[ipart], BMM_NDIM);

  dem->est.chi = v / vhc;

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    double const rho = dem->opts.part.rho / dem->est.chi;

    dem->part.m[ipart] = rho * bmm_geom_ballvol(dem->part.r[ipart], 3);
  }

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    bmm_dem_cache_j(dem, ipart);

  return true;
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
      dem->part.m[ipart] = dem->opts.part.rho * bmm_geom_ballvol(r, 3);

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
      dem->part.m[ipart] = dem->opts.part.rho * bmm_geom_ballvol(r, 3);

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
      dem->part.m[ipart] = dem->opts.part.rho * bmm_geom_ballvol(r, 3);

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
static bool bmm_dem_script_create_testblock(struct bmm_dem *const dem) {
  double const etahc = bmm_geom_ballvol(0.5, BMM_NDIM);
  double const vhc = $(bmm_prod, double)(dem->opts.box.x, BMM_NDIM);
  double const eta = dem->opts.script.params[dem->script.i].test.eta;
  double const v = vhc * (etahc / eta);

  double const vlim = vhc * eta;
  // double const rspace = bmm_ival_midpoint(dem->opts.part.rnew);
  double const rspace = dem->opts.part.rnew[1];
  double const dspace = rspace;

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
      dem->part.m[ipart] = dem->opts.part.rho * bmm_geom_ballvol(r, 3);

      for (size_t idim = 0; idim < BMM_NDIM; ++idim)
        dem->part.x[ipart][idim] = x[idim] + rspace;

      x[0] += 2.0 * rspace;

      if (x[0] + 2.0 * rspace >= dem->opts.box.x[0] + 2.0 * rspace) {
        x[0] = 0.0;
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
  double const vhc = $(bmm_prod, double)(dem->opts.box.x, BMM_NDIM);
  double const v = vhc;

  double const vlim = vhc;
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

      dem->part.role[ipart] = BMM_DEM_ROLE_FIXED;

      dem->part.r[ipart] = r;
      dem->part.m[ipart] = dem->opts.part.rho * bmm_geom_ballvol(r, 3);

      for (size_t idim = 0; idim < BMM_NDIM; ++idim)
        dem->part.x[ipart][idim] = x[idim] + rspace;

      x[0] += 2.0 * rspace / 4.0;

      if (x[0] + 2.0 * rspace >= dem->opts.box.x[0] + 2.0 * rspace) {
        break;
      }

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
      dem->part.m[ipart] = dem->opts.part.rho * bmm_geom_ballvol(r, 3);

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
      // This one is funny.
      // else dem->part.v[ipart][0] = 1.0e+2;

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
    dem->part.m[ipart] = dem->opts.part.rho * bmm_geom_ballvol(r, 3);

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

__attribute__ ((__nonnull__))
static void bmm_dem_script_perturb(struct bmm_dem *const dem) {
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    for (size_t idim = 0; idim < BMM_NDIM; ++idim)
      dem->part.x[ipart][idim] += bmm_random_get(dem->rng,
          dem->opts.part.rnew) / 2.0;
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

  double const scale = dem->opts.box.x[0];

  jpart = bmm_dem_addpart(dem);
  if (jpart == SIZE_MAX)
    return false;

  dem->part.r[jpart] = 0.015 * scale;
  dem->part.m[jpart] = dem->opts.part.rho * bmm_geom_ballvol(dem->part.r[jpart], 3);
  dem->part.x[jpart][0] += 0.45 * scale;
  dem->part.x[jpart][1] += 0.12 * scale;
  dem->part.v[jpart][0] += 0.0e+1 * scale;
  dem->part.omega[jpart] += 1.0e+4;

  jpart = bmm_dem_addpart(dem);
  if (jpart == SIZE_MAX)
    return false;

  dem->part.r[jpart] = 0.03 * scale;
  dem->part.m[jpart] = dem->opts.part.rho * bmm_geom_ballvol(dem->part.r[jpart], 3);
  dem->part.x[jpart][0] += 0.55 * scale;
  dem->part.x[jpart][1] += 0.12 * scale;
  dem->part.v[jpart][0] -= 0.0e+1 * scale;
  dem->part.omega[jpart] += 1.0e+4;

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
  dem->part.m[jpart] = dem->opts.part.rho * bmm_geom_ballvol(dem->part.r[jpart], 3);
  dem->part.x[jpart][0] += 0.45 * scale;
  dem->part.x[jpart][1] += 0.25 * scale;
  dem->part.v[jpart][0] += 3.0e+1 * scale;
  dem->part.omega[jpart] += 3.0e+4;

  jpart = bmm_dem_addpart(dem);
  if (jpart == SIZE_MAX)
    return false;

  dem->part.r[jpart] = 0.03 * scale;
  dem->part.m[jpart] = dem->opts.part.rho * bmm_geom_ballvol(dem->part.r[jpart], 3);
  dem->part.x[jpart][0] += 0.55 * scale;
  dem->part.x[jpart][1] += 0.25 * scale;
  dem->part.v[jpart][0] -= 3.0e+1 * scale;
  dem->part.omega[jpart] += 3.0e+4;

  jpart = bmm_dem_addpart(dem);
  if (jpart == SIZE_MAX)
    return false;

  dem->part.r[jpart] = 0.02 * scale;
  dem->part.m[jpart] = dem->opts.part.rho * bmm_geom_ballvol(dem->part.r[jpart], 3);
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
    case BMM_DEM_MODE_IDLE:
      break;
    case BMM_DEM_MODE_SET_DENSITY:
      bmm_dem_script_set_density(dem);

      break;
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
      // if (!bmm_dem_script_create_testpile(dem))
      //   return false;

      if (!bmm_dem_script_create_testblock(dem))
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

      dem->yield.tag = BMM_DEM_YIELD_NONE;

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
      dem->amb.tag = BMM_DEM_AMB_FAXEN;
      dem->amb.params.creeping.eta = dem->opts.script.params[dem->script.i].preset.eta;
      dem->ext.tag = BMM_DEM_EXT_NONE;
      dem->pair[BMM_DEM_CT_WEAK].cohesive = false;
      dem->pair[BMM_DEM_CT_WEAK].norm.tag = BMM_DEM_NORM_KV;
      dem->pair[BMM_DEM_CT_WEAK].norm.params.dashpot.k = dem->opts.script.params[dem->script.i].preset.kn;
      dem->pair[BMM_DEM_CT_WEAK].norm.params.dashpot.gamma = dem->opts.script.params[dem->script.i].preset.gamman;
      dem->pair[BMM_DEM_CT_WEAK].tang.tag = BMM_DEM_TANG_NONE;
      dem->pair[BMM_DEM_CT_STRONG].cohesive = true;
      dem->pair[BMM_DEM_CT_STRONG].norm.tag = BMM_DEM_NORM_KV;
      dem->pair[BMM_DEM_CT_STRONG].norm.params.dashpot.k = dem->opts.script.params[dem->script.i].preset.kn;
      dem->pair[BMM_DEM_CT_STRONG].norm.params.dashpot.gamma = dem->opts.script.params[dem->script.i].preset.gamman;
      dem->pair[BMM_DEM_CT_STRONG].tang.tag = BMM_DEM_TANG_NONE;
      dem->yield.tag = BMM_DEM_YIELD_NONE;

      break;
    case BMM_DEM_MODE_CLIP:
      bmm_dem_script_clip(dem);

      break;
    case BMM_DEM_MODE_PRESET1:
      dem->amb.tag = BMM_DEM_AMB_FAXEN;
      dem->amb.params.creeping.eta = dem->opts.script.params[dem->script.i].preset.eta2;
      dem->pair[BMM_DEM_CT_WEAK].cohesive = false;
      dem->pair[BMM_DEM_CT_WEAK].norm.tag = BMM_DEM_NORM_BSHP;
      dem->pair[BMM_DEM_CT_WEAK].norm.params.viscoel.a = dem->opts.script.params[dem->script.i].preset.a;
      dem->pair[BMM_DEM_CT_WEAK].norm.tag = BMM_DEM_NORM_KV;
      dem->pair[BMM_DEM_CT_WEAK].norm.params.dashpot.k = dem->opts.script.params[dem->script.i].preset.kn;
      dem->pair[BMM_DEM_CT_WEAK].norm.params.dashpot.gamma = dem->opts.script.params[dem->script.i].preset.gamman;
      dem->pair[BMM_DEM_CT_WEAK].tang.tag = dem->opts.script.params[dem->script.i].preset.fric;
      switch (dem->pair[BMM_DEM_CT_WEAK].tang.tag) {
        case BMM_DEM_TANG_HW:
          dem->pair[BMM_DEM_CT_WEAK].tang.params.hw.gamma = dem->opts.script.params[dem->script.i].preset.gammat;
          dem->pair[BMM_DEM_CT_WEAK].tang.params.hw.mu = dem->opts.script.params[dem->script.i].preset.mu;

          break;
        case BMM_DEM_TANG_CS:
          dem->pair[BMM_DEM_CT_WEAK].tang.params.cs.k = dem->opts.script.params[dem->script.i].preset.kt;
          dem->pair[BMM_DEM_CT_WEAK].tang.params.cs.mu = dem->opts.script.params[dem->script.i].preset.mu;

          break;
      }
      dem->pair[BMM_DEM_CT_STRONG].cohesive = true;
      dem->pair[BMM_DEM_CT_STRONG].norm.tag = BMM_DEM_NORM_KV;
      dem->pair[BMM_DEM_CT_STRONG].norm.params.dashpot.k = dem->opts.script.params[dem->script.i].preset.barkn;
      dem->pair[BMM_DEM_CT_STRONG].norm.params.dashpot.gamma = dem->opts.script.params[dem->script.i].preset.bargamman;
      dem->pair[BMM_DEM_CT_STRONG].tang.tag = BMM_DEM_TANG_BEAM;
      dem->pair[BMM_DEM_CT_STRONG].tang.params.beam.k = dem->opts.script.params[dem->script.i].preset.barkt;
      dem->pair[BMM_DEM_CT_STRONG].tang.params.beam.dk = dem->opts.script.params[dem->script.i].preset.bargammat;
      dem->yield.tag = BMM_DEM_YIELD_NONE;

      break;
    case BMM_DEM_MODE_PRESET2:
      dem->ext.tag = BMM_DEM_EXT_NONE;
      dem->amb.tag = BMM_DEM_AMB_FAXEN;
      dem->amb.params.creeping.eta = dem->opts.script.params[dem->script.i].preset.eta3;
      dem->amb.tag = BMM_DEM_AMB_NONE;
      dem->yield.tag = BMM_DEM_YIELD_ZE;
      dem->yield.params.ze.sigmacrit = dem->opts.script.params[dem->script.i].preset.sigmacrit;
      dem->yield.params.ze.taucrit = dem->opts.script.params[dem->script.i].preset.taucrit;
      dem->yield.params.ze.sigmacritt = dem->opts.script.params[dem->script.i].preset.sigmacritt;
      dem->yield.params.ze.taucritt = dem->opts.script.params[dem->script.i].preset.taucritt;

      break;
    case BMM_DEM_MODE_LINK:
      if (!bmm_dem_link(dem))
        return false;

      break;
    case BMM_DEM_MODE_FAULT:
      bmm_dem_fault(dem);

      break;
    case BMM_DEM_MODE_SEPARATE:
      // TODO Maybe implement.

      break;
    case BMM_DEM_MODE_EXPORT:
      if (dem->opts.script.params[dem->script.i].expr.entropic) {
        dump_raddist_etc(dem);

        if (!export_s(dem) || !export_chi(dem)) {
          BMM_TLE_EXTS(BMM_TLE_NUM_UNKNOWN, "Nope");

          return false;
        }
      }

      if (!export_x(dem) || !export_phi(dem) ||
          !export_r(dem) || !export_c(dem) ||
          !export_f(dem) || !export_p(dem)) {
        BMM_TLE_EXTS(BMM_TLE_NUM_UNKNOWN, "Big nope");

        return false;
      }

      break;
    case BMM_DEM_MODE_PRECRUNCH:
      for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
        if (dem->part.x[ipart][1] > dem->opts.box.x[1] -
            dem->opts.script.params[dem->script.i].precrunch.nlayer *
            2.0 * bmm_ival_midpoint(dem->opts.part.rnew))
          dem->part.role[ipart] = BMM_DEM_ROLE_DRIVEN;

      for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
        if (dem->part.x[ipart][1] <
            dem->opts.script.params[dem->script.i].precrunch.nlayer *
            2.0 * bmm_ival_midpoint(dem->opts.part.rnew))
          dem->part.role[ipart] = BMM_DEM_ROLE_FIXED;

      break;
    case BMM_DEM_MODE_ZEROEST:
      dem->est.eambdis = 0.0;
      dem->est.edrivtang = 0.0;
      dem->est.edrivnorm = 0.0;
      dem->est.ebond = 0.0;
      dem->est.eyieldis = 0.0;
      dem->est.ewcontdis = 0.0;
      dem->est.escontdis = 0.0;
      dem->est.hwgamma = 0;
      dem->est.hwmu = 0;
      dem->est.csk = 0;
      dem->est.csmu = 0;

      break;
    case BMM_DEM_MODE_CRUNCH:
      {
        static size_t stage = 0;
        if (stage != dem->script.i) {
          if (stage == 0)
            for (size_t idim = 0; idim < BMM_NDIM; ++idim)
              dem->script.state.crunch.fdrive[idim] = 0.0;

          stage = dem->script.i;

          dem->ext.tag = BMM_DEM_EXT_DRIVE;

          dem->script.state.crunch.ndrive = 0;
          for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
            if (dem->part.role[ipart] == BMM_DEM_ROLE_DRIVEN)
              ++dem->script.state.crunch.ndrive;

          dem->script.state.crunch.nfix = 0;
          for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
            if (dem->part.role[ipart] == BMM_DEM_ROLE_FIXED) {
              ++dem->script.state.crunch.nfix;

              for (size_t idim = 0; idim < BMM_NDIM; ++idim)
                dem->part.a[ipart][idim] = 0.0;

              dem->part.alpha[ipart] = 0.0;

              for (size_t idim = 0; idim < BMM_NDIM; ++idim)
                dem->part.v[ipart][idim] = 0.0;

              dem->part.omega[ipart] = 0.0;
            }
        }
      }

      break;
    case BMM_DEM_MODE_STORE:
      {
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
      }

      break;
    case BMM_DEM_MODE_LOAD:
      {
        FILE *const stream = fopen("a.out", "r");
        if (stream == NULL) {
          BMM_TLE_STDS();
          abort();
        }
        if (fread(dem, sizeof *dem, 1, stream) != 1) {
          BMM_TLE_STDS();
          abort();
        }
        if (fclose(stream) != 0) {
          BMM_TLE_STDS();
          abort();
        }

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
      }

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
  stream = fopen("est.data", "w");
  if (stream == NULL) {
    BMM_TLE_STDS();

    return false;
  }

  return true;
}

static bool garbage(struct bmm_dem const *const dem) {
  double const eambdis = dem->est.eambdis;
  double const epotext = dem->est.epotext_d;
  double const eklin = dem->est.eklin_d;
  double const ekrot = dem->est.ekrot_d;
  double const ewcont = dem->est.ewcont;
  double const escont = dem->est.escont;
  double const edrivnorm = dem->est.edrivnorm;
  double const edrivtang = dem->est.edrivtang;
  double const ebond = dem->est.ebond;
  double const eyieldis = dem->est.eyieldis;
  double const ewcontdis = dem->est.ewcontdis;
  double const escontdis = dem->est.escontdis;
  double const pos = eambdis + epotext + eklin + ekrot + ewcont + escont;
  double const dis = eyieldis + ewcontdis + escontdis;
  double const neg = edrivnorm + edrivtang + ebond;
  double const eee = pos + dis - neg;

  if (dem->opts.script.mode[dem->script.i] == BMM_DEM_MODE_CRUNCH &&
      dem->opts.script.params[dem->script.i].crunch.measure)
    if (fprintf(stream, "%g %g %g %g %g %g %g %g %g %g\n",
          dem->time.t,
          pos, dis, neg,
          dem->est.mueff, dem->est.mueffb,
          dem->est.vdriv[0], dem->est.vdriv[1],
          dem->script.state.crunch.fdrive[0], dem->script.state.crunch.fdrive[1]) < 0 ||
        fflush(stream) == EOF) {
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

    if (!bmm_dem_puts(dem, BMM_MSG_NUM_EST))
      return false;

    if (!garbage(dem)) {
      BMM_TLE_EXTS(BMM_TLE_NUM_UNKNOWN, "Nope");

      abort();
    }
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

    // TODO This is lost here.

    if (fprintf(stderr, "P-Wave Speed: %g\n", bmm_dem_pwavec(dem)) < 0) {
      BMM_TLE_STDS();

      return false;
    }

    if (fprintf(stderr, "S-Wave Speed: %g\n", bmm_dem_swavec(dem)) < 0) {
      BMM_TLE_STDS();

      return false;
    }

    switch (dem->pair[BMM_DEM_CT_WEAK].tang.tag) {
      case BMM_DEM_TANG_HW:
        if (fprintf(stderr, "HW Dynamic Fraction: %g\n",
              (double) dem->est.hwmu /
              (double) (dem->est.hwgamma + dem->est.hwmu)) < 0) {
          BMM_TLE_STDS();

          return false;
        }

        break;
      case BMM_DEM_TANG_CS:
        if (fprintf(stderr, "CS Dynamic Fraction: %g\n",
              (double) dem->est.csmu /
              (double) (dem->est.csk + dem->est.csmu)) < 0) {
          BMM_TLE_STDS();

          return false;
        }

        break;
    }
  }

  return true;
}

static bool bmm_dem_run_(struct bmm_dem *const dem) {
  int const sigs[] = {SIGUSR1, SIGUSR2, SIGINT, SIGQUIT, SIGTERM, SIGPIPE};
  if (bmm_sig_register(sigs, nmembof(sigs)) != SIZE_MAX) {
    BMM_TLE_STDS();

    return false;
  }

  if (!pregarbage(dem))
    BMM_TLE_EXTS(BMM_TLE_NUM_UNKNOWN, "Nope");

  for ever {
    int signum;
    if (bmm_sig_use(&signum))
      switch (signum) {
        case SIGUSR1:
        case SIGUSR2:
          fprintf(stderr, "Time: %g, Script: %zu / %zu\n",
              dem->time.t, dem->script.i + 1, dem->opts.script.n);

          break;
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

  bmm_dem_trap_off(dem);

  if (!postgarbage(dem))
    BMM_TLE_EXTS(BMM_TLE_NUM_UNKNOWN, "Nope");

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
