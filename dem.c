#include <gsl/gsl_rng.h>
#include <math.h>
#include <signal.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "conf.h"
#include "cpp.h"
#include "dem.h"
#include "fp.h"
#include "geom.h"
#include "geom2d.h"
#include "io.h"
#include "neigh.h"
#include "msg.h"
#include "random.h"
#include "sig.h"
#include "size.h"
#include "tle.h"

size_t bmm_dem_script_addstage(struct bmm_dem_opts* const opts) {
  size_t const istage = opts->script.n;

  if (istage >= BMM_MSTAGE)
    return SIZE_MAX;

  ++opts->script.n;

  opts->script.mode[istage] = BMM_DEM_MODE_IDLE;
  opts->script.tspan[istage] = 0.0;
  opts->script.dt[istage] = 0.0;

  return istage;
}

bool bmm_dem_script_ongoing(struct bmm_dem const* const dem) {
  return dem->script.i < dem->opts.script.n;
}

bool bmm_dem_script_trans(struct bmm_dem* const dem) {
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
static void bmm_dem_cache_j(struct bmm_dem* const dem,
    size_t const ipart) {
  dem->cache.j[ipart] = dem->part.jred[ipart] *
    dem->part.m[ipart] * bmm_fp_sq(dem->part.r[ipart]);
}

/// The call `bmm_dem_cache_x(dem, ipart)`
/// caches the position of the particle `ipart`
/// in the simulation `dem`.
__attribute__ ((__nonnull__))
static void bmm_dem_cache_x(struct bmm_dem* const dem,
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
static void bmm_dem_cache_ijcell(struct bmm_dem* const dem,
    size_t const ipart) {
  for (size_t idim = 0; idim < BMM_NDIM; ++idim) {
    size_t const i = dem->opts.cache.ncell[idim] - 1;

    double const a = 1.0;
    double const b = (double) i;
    double const t = bmm_fp_lerp(dem->cache.x[ipart][idim],
        0.0, dem->opts.box.x[idim], a, b);

    if (t < a)
      dem->cache.ijcell[ipart][idim] = 0;
    else if (t >= b)
      dem->cache.ijcell[ipart][idim] = i;
    else {
      size_t const j = (size_t) t;

      // The equivalent for a signed type would be
      // `dynamic_assert(j >= 0 && j < i, "Invalid truncation")`.
      dynamic_assert(j < i, "Invalid truncation");

      dem->cache.ijcell[ipart][idim] = j;
    }
  }
}

/// The call `bmm_dem_cache_icell(dem, ipart)`
/// caches the neighbor cell index of the particle `ipart`
/// in the simulation `dem`.
/// The neighbor cell index vectors need to be cached first
/// by calling `bmm_dem_cache_ijcell`.
__attribute__ ((__nonnull__))
static void bmm_dem_cache_icell(struct bmm_dem* const dem,
    size_t const ipart) {
  dem->cache.icell[ipart] = bmm_size_unhcd(dem->cache.ijcell[ipart],
      BMM_NDIM, dem->opts.cache.ncell);
}

/// The call `bmm_dem_cache_clrparts(dem)`
/// clears the neighbor cell index mapping cache
/// in the simulation `dem`.
__attribute__ ((__nonnull__))
static void bmm_dem_cache_clrparts(struct bmm_dem* const dem) {
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
static bool bmm_dem_cache_addpart(struct bmm_dem* const dem,
    size_t const ipart) {
  size_t const icell = dem->cache.icell[ipart];

  if (dem->cache.part[icell].n >= nmembof(dem->cache.part[icell].i)) {
    BMM_TLE_EXTS(BMM_TLE_NUM_UNKNOWN, "Neighbor cell capacity exhausted");

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
static bool bmm_dem_cache_eligible(struct bmm_dem const* const dem,
    size_t const ipart, size_t const jpart) {
  if (dem->cache.icell[ipart] == dem->cache.icell[jpart] && jpart <= ipart)
    return false;

  if (bmm_geom2d_cpdist2(dem->cache.x[ipart], dem->cache.x[jpart],
        dem->opts.box.x, dem->opts.box.per) >
      bmm_fp_sq(dem->opts.cache.rcutoff))
    return false;

  return true;
}

/// The call `bmm_dem_cache_clrneighs(dem)`
/// clears the neighbor cache
/// in the simulation `dem`.
__attribute__ ((__nonnull__))
static void bmm_dem_cache_clrneighs(struct bmm_dem* const dem) {
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
static bool bmm_dem_cache_addfrom(struct bmm_dem* const dem,
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
          BMM_TLE_EXTS(BMM_TLE_NUM_UNKNOWN, "Neighbor capacity exhausted");

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
static bool bmm_dem_cache_addto(struct bmm_dem* const dem,
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
          BMM_TLE_EXTS(BMM_TLE_NUM_UNKNOWN, "Neighbor capacity exhausted");

          return false;
        }

        dem->cache.neigh[jpart].i[dem->cache.neigh[jpart].n] = ipart;
        ++dem->cache.neigh[jpart].n;
      }
    }
  }

  return true;
}

bool bmm_dem_cache_build(struct bmm_dem* const dem) {
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

size_t bmm_dem_addpart(struct bmm_dem* const dem) {
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

  for (size_t istep = 0; istep < nmembof(dem->part.aprev); ++istep)
    for (size_t idim = 0; idim < BMM_NDIM; ++idim)
      dem->part.aprev[istep][ipart][idim] = 0.0;

  dem->part.phi[ipart] = 0.0;
  dem->part.omega[ipart] = 0.0;
  dem->part.alpha[ipart] = 0.0;

  for (size_t istep = 0; istep < nmembof(dem->part.alphaprev); ++istep)
    dem->part.alphaprev[istep][ipart] = 0.0;

  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    dem->part.f[ipart][idim] = 0.0;

  dem->part.tau[ipart] = 0.0;

  dem->cache.stale = true;

  return ipart;
}

/// The call `bmm_dem_reassign(dem, ipart, jpart)`
/// reassigns the particle `jpart` to `ipart`.
__attribute__ ((__nonnull__))
static void bmm_dem_reassign(struct bmm_dem* const dem,
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

  for (size_t istep = 0; istep < nmembof(dem->part.aprev); ++istep)
    for (size_t idim = 0; idim < BMM_NDIM; ++idim)
      dem->part.aprev[istep][ipart][idim] =
        dem->part.aprev[istep][jpart][idim];

  dem->part.phi[ipart] = dem->part.phi[jpart];
  dem->part.omega[ipart] = dem->part.omega[jpart];
  dem->part.alpha[ipart] = dem->part.alpha[jpart];

  for (size_t istep = 0; istep < nmembof(dem->part.alphaprev); ++istep)
    dem->part.alphaprev[istep][ipart] = dem->part.alphaprev[istep][jpart];

  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    dem->part.f[ipart][idim] = dem->part.f[jpart][idim];

  dem->part.tau[ipart] = dem->part.tau[jpart];

  dem->link.part[ipart].n = dem->link.part[jpart].n;

  for (size_t ilink = 0; ilink < dem->link.part[jpart].n; ++ilink)
    dem->link.part[ipart].i[ilink] = dem->link.part[jpart].i[ilink];

  for (size_t ilink = 0; ilink < dem->link.part[jpart].n; ++ilink)
    dem->link.part[ipart].rrest[ilink] = dem->link.part[jpart].rrest[ilink];

  for (size_t ilink = 0; ilink < dem->link.part[jpart].n; ++ilink)
    for (size_t iend = 0; iend < nmembof(dem->link.part[ipart].phirest[ilink]); ++iend)
      dem->link.part[ipart].phirest[ilink][iend] =
        dem->link.part[jpart].phirest[ilink][iend];

  for (size_t ilink = 0; ilink < dem->link.part[jpart].n; ++ilink)
    dem->link.part[ipart].rlim[ilink] = dem->link.part[jpart].rlim[ilink];

  for (size_t ilink = 0; ilink < dem->link.part[jpart].n; ++ilink)
    dem->link.part[ipart].philim[ilink] = dem->link.part[jpart].philim[ilink];
}

void bmm_dem_rempart(struct bmm_dem* const dem,
    size_t const ipart) {
  --dem->part.n;

  size_t const jpart = dem->part.n;

  if (jpart != ipart)
    bmm_dem_reassign(dem, ipart, jpart);

  dem->cache.stale = true;
}

void bmm_dem_force_pair(struct bmm_dem* const dem,
    size_t const ipart, size_t const jpart) {
  double xdiff[BMM_NDIM];
  bmm_geom2d_cpdiff(xdiff, dem->part.x[ipart], dem->part.x[jpart],
      dem->opts.box.x, dem->opts.box.per);

  double const d2 = bmm_geom2d_norm2(xdiff);

  double const r = dem->part.r[ipart] + dem->part.r[jpart];
  double const r2 = bmm_fp_sq(r);

  // TODO Investigate.

  // TODO !!
  if (d2 > r2 || d2 == 0.0)
    return;

  double const d = sqrt(d2);

  double xnorm[BMM_NDIM];
  bmm_geom2d_scale(xnorm, xdiff, 1.0 / d);

  double xtang[BMM_NDIM];
  bmm_geom2d_rperp(xtang, xnorm);

  double nvdiff[BMM_NDIM];
  bmm_geom2d_diff(nvdiff, dem->part.v[ipart], dem->part.v[jpart]);
  // TODO ??
  double vdiff[BMM_NDIM];
  bmm_geom2d_scale(vdiff, nvdiff, -1.0);

  double const xi = r - d;
  double const dotxi = bmm_geom2d_dot(vdiff, xnorm);
  double const vtang = bmm_geom2d_dot(vdiff, xtang) +
    dem->part.r[ipart] * dem->part.omega[ipart] +
    dem->part.r[jpart] * dem->part.omega[jpart];

  double fn = 0.0;
  switch (dem->fnorm) {
    case BMM_DEM_FNORM_DASHPOT:
      fn = fmax(0.0, dem->opts.part.y * xi +
          dem->norm.params.dashpot.gamma * dotxi);

      break;
  }

  double fdiff[2];
  bmm_geom2d_scale(fdiff, xnorm, -fn);

  bmm_geom2d_addto(dem->part.f[ipart], fdiff);

  double fdiff2[2];
  bmm_geom2d_scale(fdiff2, fdiff, -1.0);

  bmm_geom2d_addto(dem->part.f[jpart], fdiff2);

  double ft = 0.0;
  switch (dem->ftang) {
    case BMM_DEM_FTANG_HW:
      ft = -copysign(fmin(dem->tang.params.hw.gamma * fabs(vtang),
            dem->tang.params.hw.mu * fn), vtang);

      break;
  }

  bmm_geom2d_scale(fdiff, xtang, ft);

  bmm_geom2d_addto(dem->part.f[ipart], fdiff);

  dem->part.tau[ipart] += ft * dem->part.r[ipart];

  bmm_geom2d_scale(fdiff2, fdiff, -1.0);

  bmm_geom2d_addto(dem->part.f[jpart], fdiff2);

  dem->part.tau[jpart] += ft * dem->part.r[jpart];
}

void bmm_dem_force_creeping(struct bmm_dem* const dem,
    size_t const ipart) {
  double const v = bmm_geom2d_norm(dem->part.v[ipart]);

  if (v == 0.0)
    return;

  double vunit[BMM_NDIM];
  bmm_geom2d_scale(vunit, dem->part.v[ipart], 1.0 / v);

  double const f = -3.0 * M_2PI * dem->ambient.params.creeping.mu *
    dem->part.r[ipart] * v;

  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    dem->part.f[ipart][idim] += f * vunit[idim];

  double const tau = -4.0 * M_2PI * dem->ambient.params.creeping.mu *
    bmm_fp_cb(dem->part.r[ipart]);

  dem->part.tau[ipart] += tau * dem->part.omega[ipart];
}

// TODO Flatten these to allow loop-invariant code motion.
void bmm_dem_force_ambient(struct bmm_dem* const dem, size_t const ipart) {
  switch (dem->famb) {
    case BMM_DEM_FAMB_CREEPING:
      bmm_dem_force_creeping(dem, ipart);

      break;
    case BMM_DEM_FAMB_QUAD:

      break;
    case BMM_DEM_FAMB_CORR:

      break;
  }
}

void bmm_dem_force_link(struct bmm_dem* const dem,
    size_t const ipart, size_t const jpart, size_t const ilink) {
  // TODO This.
  return;
}

void bmm_dem_force_external(struct bmm_dem* const dem, size_t const ipart) {
  switch (dem->fext) {
    case BMM_DEM_FEXT_ABS:
      dem->part.f[ipart][1] += copysign(dem->ext.params.abs.fcohes,
          dem->opts.box.x[1] / 2.0 - dem->part.x[ipart][1]);

      break;
    case BMM_DEM_FEXT_HARM:
      dem->part.f[ipart][1] += dem->ext.params.harm.kcohes *
        (dem->opts.box.x[1] / 2.0 - dem->part.x[ipart][1]);

      break;
  }
}

void bmm_dem_force(struct bmm_dem* const dem) {
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    for (size_t idim = 0; idim < BMM_NDIM; ++idim)
      dem->part.f[ipart][idim] = 0.0;

    dem->part.tau[ipart] = 0.0;
  }

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    bmm_dem_force_ambient(dem, ipart);

  switch (dem->caching) {
    case BMM_DEM_CACHING_NONE:
      for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
        for (size_t jpart = ipart + 1; jpart < dem->part.n; ++jpart)
          bmm_dem_force_pair(dem, ipart, jpart);

      break;
    case BMM_DEM_CACHING_NEIGH:
      for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
        for (size_t ineigh = 0; ineigh < dem->cache.neigh[ipart].n; ++ineigh)
          bmm_dem_force_pair(dem, ipart, dem->cache.neigh[ipart].i[ineigh]);

      break;
  }

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    for (size_t ilink = 0; ilink < dem->link.part[ipart].n; ++ilink)
      bmm_dem_force_link(dem, ipart, dem->link.part[ipart].i[ilink], ilink);

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    bmm_dem_force_external(dem, ipart);

  // TODO Calculate force feedback from the residuals.
}

void bmm_dem_integ_euler(struct bmm_dem* const dem) {
  double const dt = dem->opts.script.dt[dem->script.i];

  if (dt == 0.0)
    return;

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    double const m = dem->part.m[ipart];
    double const j = dem->cache.j[ipart];

    for (size_t idim = 0; idim < BMM_NDIM; ++idim)
      dem->part.a[ipart][idim] = dem->part.f[ipart][idim] / m;

    dem->part.alpha[ipart] = dem->part.tau[ipart] / j;
  }

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    for (size_t idim = 0; idim < BMM_NDIM; ++idim) {
      dem->part.x[ipart][idim] += dem->part.v[ipart][idim] * dt;
      dem->part.v[ipart][idim] += dem->part.a[ipart][idim] * dt;

      if (dem->opts.box.per[idim])
        dem->part.x[ipart][idim] = bmm_fp_uwrap(dem->part.x[ipart][idim],
            dem->opts.box.x[idim]);
    }

    dem->part.phi[ipart] += dem->part.omega[ipart] * dt;
    dem->part.omega[ipart] += dem->part.alpha[ipart] * dt;
  }
}

void bmm_dem_integ_taylor3(struct bmm_dem* const dem) {
  double const dt = dem->opts.script.dt[dem->script.i];

  if (dt == 0.0)
    return;

  double const dt2 = dt * dt;
  double const dt3 = dt2 * dt;

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    for (size_t idim = 0; idim < BMM_NDIM; ++idim)
      dem->part.aprev[0][ipart][idim] = dem->part.a[ipart][idim];

    dem->part.alphaprev[0][ipart] = dem->part.alpha[ipart];
  }

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    double const m = dem->part.m[ipart];
    double const j = dem->cache.j[ipart];

    for (size_t idim = 0; idim < BMM_NDIM; ++idim)
      dem->part.a[ipart][idim] = dem->part.f[ipart][idim] / m;

    dem->part.alpha[ipart] = dem->part.tau[ipart] / j;
  }

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    for (size_t idim = 0; idim < BMM_NDIM; ++idim) {
      double const b = (dem->part.a[ipart][idim] -
          dem->part.aprev[0][ipart][idim]) / dt;

      dem->part.x[ipart][idim] +=
        (1.0 / bmm_fp_fact(1)) * dem->part.v[ipart][idim] * dt +
        (1.0 / bmm_fp_fact(2)) * dem->part.a[ipart][idim] * dt2 +
        (1.0 / bmm_fp_fact(3)) * b * dt3;
      dem->part.v[ipart][idim] +=
        (1.0 / bmm_fp_fact(1)) * dem->part.a[ipart][idim] * dt +
        (1.0 / bmm_fp_fact(2)) * b * dt2;

      if (dem->opts.box.per[idim])
        dem->part.x[ipart][idim] = bmm_fp_uwrap(dem->part.x[ipart][idim],
            dem->opts.box.x[idim]);
    }

    double const beta = (dem->part.alpha[ipart] -
        dem->part.alphaprev[0][ipart]) / dt;

    dem->part.phi[ipart] +=
      (1.0 / bmm_fp_fact(1)) * dem->part.omega[ipart] * dt;
      (1.0 / bmm_fp_fact(2)) * dem->part.alpha[ipart] * dt2 +
      (1.0 / bmm_fp_fact(3)) * beta * dt3;
    dem->part.omega[ipart] +=
      (1.0 / bmm_fp_fact(1)) * dem->part.alpha[ipart] * dt +
      (1.0 / bmm_fp_fact(2)) * beta * dt2;
  }
}

void bmm_dem_integ_gear4pred(struct bmm_dem* const dem) {
  double const dt = dem->opts.script.dt[dem->script.i];

  if (dt == 0.0)
    return;

  double const dt2 = dt * dt;
  double const dt3 = dt2 * dt;

  // TODO Same as part of `taylor3`.
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    for (size_t idim = 0; idim < BMM_NDIM; ++idim)
      dem->pred.gear.aprev[0][ipart][idim] = dem->part.a[ipart][idim];

    dem->pred.gear.alphaprev[0][ipart] = dem->part.alpha[ipart];
  }

  // TODO Almost the same as part of `taylor3`.
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    for (size_t idim = 0; idim < BMM_NDIM; ++idim) {
      // double const b = (dem->part.a[ipart][idim] - dem->part.aprev[0][ipart][idim]) / dt;
      double const b = dem->pred.gear.b[ipart][idim];

      dem->part.x[ipart][idim] += dem->part.v[ipart][idim] * dt +
        (1.0 / 2.0) * dem->part.a[ipart][idim] * dt2 +
        (1.0 / 6.0) * b * dt3;
      dem->part.v[ipart][idim] += dem->part.a[ipart][idim] * dt +
        (1.0 / 2.0) * b * dt2;
      dem->part.a[ipart][idim] += b * dt;
      dem->pred.gear.b[ipart][idim] = b;

      if (dem->opts.box.per[idim])
        dem->part.x[ipart][idim] = bmm_fp_uwrap(dem->part.x[ipart][idim],
            dem->opts.box.x[idim]);
    }

    // double const beta = (dem->part.alpha[ipart] - dem->part.alphaprev[0][ipart]) / dt;
    double const beta = dem->pred.gear.beta[ipart];

    dem->part.phi[ipart] += dem->part.omega[ipart] * dt;
      (1.0 / 2.0) * dem->part.alpha[ipart] * dt2 +
      (1.0 / 6.0) * beta * dt3;
    dem->part.omega[ipart] += dem->part.alpha[ipart] * dt +
      (1.0 / 2.0) * beta * dt2;
    dem->part.alpha[ipart] += beta * dt;
    dem->pred.gear.beta[ipart] = beta;
  }
}

void bmm_dem_integ_gear4corr(struct bmm_dem* const dem) {
  static double const c[] = {1.0 / 6.0, 5.0 / 6.0, 1.0 / 1.0, 1.0 / 3.0};

  double const dt = dem->opts.script.dt[dem->script.i];

  if (dt == 0.0)
    return;

  double const dt2 = dt * dt;
  double const dt3 = dt2 * dt;

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    for (size_t idim = 0; idim < BMM_NDIM; ++idim)
      dem->pred.gear.a[ipart][idim] = dem->part.a[ipart][idim];

    dem->pred.gear.alpha[ipart] = dem->part.alpha[ipart];
  }

  // TODO Same as part of `taylor3`.
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    double const m = dem->part.m[ipart];
    double const j = dem->cache.j[ipart];

    for (size_t idim = 0; idim < BMM_NDIM; ++idim)
      dem->part.a[ipart][idim] = dem->part.f[ipart][idim] / m;

    dem->part.alpha[ipart] = dem->part.tau[ipart] / j;
  }

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    for (size_t idim = 0; idim < BMM_NDIM; ++idim) {
      double const da = dem->part.a[ipart][idim] -
          dem->pred.gear.a[ipart][idim];

      dem->part.x[ipart][idim] +=
        bmm_fp_fact(0) * (1.0 / 2.0) * c[0] * da * dt2;
      dem->part.v[ipart][idim] +=
        bmm_fp_fact(1) * (1.0 / 2.0) * c[1] * da * dt;
      dem->part.a[ipart][idim] +=
        bmm_fp_fact(2) * (1.0 / 2.0) * c[2] * da;
      dem->pred.gear.b[ipart][idim] +=
        bmm_fp_fact(3) * (1.0 / 2.0) * c[3] * da / dt;

      if (dem->opts.box.per[idim])
        dem->part.x[ipart][idim] = bmm_fp_uwrap(dem->part.x[ipart][idim],
            dem->opts.box.x[idim]);
    }

    double const dalpha = dem->part.alpha[ipart] - dem->pred.gear.alpha[ipart];

    dem->part.phi[ipart] +=
      bmm_fp_fact(0) * (1.0 / 2.0) * c[0] * dalpha * dt2;
    dem->part.omega[ipart] +=
      bmm_fp_fact(1) * (1.0 / 2.0) * c[1] * dalpha * dt;
    dem->part.alpha[ipart] +=
      bmm_fp_fact(2) * (1.0 / 2.0) * c[2] * dalpha;
    dem->pred.gear.beta[ipart] +=
      bmm_fp_fact(3) * (1.0 / 2.0) * c[3] * dalpha / dt;
  }

  // TODO Same as part of `taylor3`.
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    for (size_t idim = 0; idim < BMM_NDIM; ++idim)
      dem->part.aprev[0][ipart][idim] = dem->pred.gear.aprev[0][ipart][idim];

    dem->part.alphaprev[0][ipart] = dem->pred.gear.alphaprev[0][ipart];
  }
}

void bmm_dem_stab(struct bmm_dem* const dem) {
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    dem->part.phi[ipart] = bmm_fp_uwrap(dem->part.phi[ipart], M_2PI);
}

void bmm_dem_predict(struct bmm_dem* const dem) {
  switch (dem->integ) {
    case BMM_DEM_INTEG_GEAR:
      bmm_dem_integ_gear4pred(dem);

      break;
  }
}

void bmm_dem_correct(struct bmm_dem* const dem) {
  switch (dem->integ) {
    case BMM_DEM_INTEG_EULER:
      bmm_dem_integ_euler(dem);

      break;
    case BMM_DEM_INTEG_TAYLOR3:
      bmm_dem_integ_taylor3(dem);

      break;
    case BMM_DEM_INTEG_GEAR:
      bmm_dem_integ_gear4corr(dem);

      break;
  }
}

bool bmm_dem_link_pair(struct bmm_dem* const dem,
    size_t const ipart, size_t const jpart) {
  double xdiff[BMM_NDIM];
  bmm_geom2d_cpdiff(xdiff, dem->part.x[ipart], dem->part.x[jpart],
      dem->opts.box.x, dem->opts.box.per);

  double const d2 = bmm_geom2d_norm2(xdiff);

  double const r = dem->part.r[ipart] + dem->part.r[jpart];
  double const r2 = bmm_fp_sq(r);

  if (d2 > r2 * dem->opts.link.ccrlink)
    return false;

  if (dem->link.part[ipart].n >= nmembof(dem->link.part[ipart].i))
    return false;

  double const d = sqrt(d2);

  double const rrest = d * dem->opts.link.cshlink;

  // TODO Check sign.
  switch (dem->flink) {
    case BMM_DEM_FLINK_BEAM:
      {
        double const phi = bmm_geom2d_dir(xdiff);

        dem->link.part[ipart].phirest[dem->link.part[ipart].n][0] =
          dem->part.phi[ipart] - phi;

        dem->link.part[ipart].phirest[dem->link.part[ipart].n][1] =
          dem->part.phi[jpart] - bmm_geom2d_redir(phi);

        double const crlim = gsl_rng_uniform(dem->rng) *
          (dem->opts.link.crlim[1] - dem->opts.link.crlim[0]) +
          dem->opts.link.crlim[0];

        dem->link.part[ipart].rlim[dem->link.part[ipart].n] =
          crlim * rrest;

        double const cphilim = gsl_rng_uniform(dem->rng) *
          (dem->opts.link.cphilim[1] - dem->opts.link.cphilim[0]) +
          dem->opts.link.cphilim[0];

        dem->link.part[ipart].philim[dem->link.part[ipart].n] =
          cphilim * M_2PI;
      }

      break;
  }

  dem->link.part[ipart].rrest[dem->link.part[ipart].n] = rrest;
  dem->link.part[ipart].i[dem->link.part[ipart].n] = jpart;
  ++dem->link.part[ipart].n;

  return true;
}

// TODO Check triangulation quality and compare with Delaunay.
bool bmm_dem_link(struct bmm_dem* const dem) {
  switch (dem->caching) {
    case BMM_DEM_CACHING_NONE:
      for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
        for (size_t jpart = ipart + 1; jpart < dem->part.n; ++jpart)
          bmm_dem_link_pair(dem, ipart, jpart);

      break;
    case BMM_DEM_CACHING_NEIGH:
      for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
        for (size_t ineigh = 0; ineigh < dem->cache.neigh[ipart].n; ++ineigh)
          bmm_dem_link_pair(dem, ipart, dem->cache.neigh[ipart].i[ineigh]);

      break;
  }

  return true;
}

bool bmm_dem_unlink(struct bmm_dem* const dem) {
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    dem->link.part[ipart].n = 0;

  return true;
}

/// The call `bmm_dem_est_ekin(dem)`
/// returns the total kinetic energy of the particles
/// in the simulation `dem`.
__attribute__ ((__nonnull__, __pure__))
double bmm_dem_est_ekin(struct bmm_dem const* const dem) {
  double e = 0.0;

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    for (size_t idim = 0; idim < 2; ++idim)
      e += dem->part.m[ipart] * bmm_fp_sq(dem->part.v[ipart][idim]);

    e += dem->cache.j[ipart] * bmm_fp_sq(dem->part.omega[ipart]);
  }

  return (1.0 / 2.0) * e;
}

/// The call `bmm_dem_est_mass(dem)`
/// returns the total mass of the particles
/// in the simulation `dem`.
__attribute__ ((__nonnull__, __pure__))
double bmm_dem_est_mass(struct bmm_dem* const dem) {
  double m = 0.0;

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    m += dem->part.m[ipart];

  return m;
}

/// The call `bmm_dem_est_center(pxcenter, dem)`
/// sets `pxcenter` to the center of the bounding box
/// in the simulation `dem`.
__attribute__ ((__nonnull__))
void bmm_dem_est_center(double* const pxcenter,
    struct bmm_dem* const dem) {
  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    pxcenter[idim] = dem->opts.box.x[idim] / 2.0;
}

/// The call `bmm_dem_est_com(pxcom, dem)`
/// sets `pxcom` to the center of mass of the particles
/// in the simulation `dem`.
__attribute__ ((__nonnull__))
void bmm_dem_est_com(double* const pxcom,
    struct bmm_dem* const dem) {
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
double bmm_dem_est_cor(struct bmm_dem const* const dem) {
  double e = 0.0;

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    double const mred =
      dem->part.m[ipart] * dem->part.m[ipart] /
      (dem->part.m[ipart] + dem->part.m[ipart]);
    e += exp(-M_PI * dem->norm.params.dashpot.gamma / (2.0 * mred) /
        sqrt(dem->opts.part.y / mred -
          bmm_fp_sq(dem->norm.params.dashpot.gamma / (2.0 * mred))));
  }

  return e / (double) dem->part.n;
}

void bmm_dem_opts_def(struct bmm_dem_opts* const opts) {
  // This is here just to help Valgrind and cover up my mistakes.
  (void) memset(opts, 0, sizeof *opts);

  opts->verbose = false;

  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    opts->box.x[idim] = 1.0;

  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    opts->box.per[idim] = false;

  opts->time.istab = 1000;

  opts->part.y = 1.0;
  opts->part.rnew[0] = 1.0;
  opts->part.rnew[1] = 1.0;

  opts->link.ccrlink = 1.2;
  opts->link.cshlink = 0.8;
  opts->link.ktens = 1.0;
  opts->link.kshear = 1.0;
  opts->link.crlim[0] = 1.0;
  opts->link.crlim[1] = 1.0;
  opts->link.cphilim[0] = 1.0;
  opts->link.cphilim[1] = 1.0;

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

  opts->cache.rcutoff = (double) INFINITY;

  for (size_t idim = 0; idim < BMM_NDIM; ++idim) {
    opts->cache.rcutoff = fmin(opts->cache.rcutoff,
        opts->box.x[idim] / (double) (opts->cache.ncell[idim] - 2));

    // TODO Use this somewhere.
    /*
    dynamic_assert(opts->cache.rcutoff <= opts->box.x[idim] /
        (double) ((opts->cache.ncell[idim] - 2) * 2),
        "Neighbor cells too small");
    */
  }
}

void bmm_dem_def(struct bmm_dem* const dem,
    struct bmm_dem_opts const* const opts) {
  // This is here just to help Valgrind and cover up my mistakes.
  (void) memset(dem, 0, sizeof *dem);

  dem->opts = *opts;

  dem->integ = BMM_DEM_INTEG_EULER;
  dem->integ = BMM_DEM_INTEG_TAYLOR3;
  dem->integ = BMM_DEM_INTEG_GEAR;
  dem->caching = BMM_DEM_CACHING_NONE;
  dem->caching = BMM_DEM_CACHING_NEIGH;
  dem->fext = BMM_DEM_FEXT_NONE;
  dem->famb = BMM_DEM_FAMB_CREEPING;
  dem->fnorm = BMM_DEM_FNORM_DASHPOT;
  dem->ftang = BMM_DEM_FTANG_HW;
  dem->flink = BMM_DEM_FLINK_BEAM;

  dem->ambient.params.creeping.mu = 1.0;
  dem->norm.params.dashpot.gamma = 1.0e+1;
  dem->tang.params.hw.gamma = 1.0;
  dem->tang.params.hw.mu = 1.0;

  dem->time.t = 0.0;
  dem->time.istep = 0;

  dem->part.n = 0;
  dem->part.lnew = 0;

  for (size_t ipart = 0; ipart < nmembof(dem->link.part); ++ipart)
    dem->link.part[ipart].n = 0;

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

static bool msg_write(void const* buf, size_t const n,
    __attribute__ ((__unused__)) void* const ptr) {
  return bmm_io_writeout(buf, n);
}

size_t bmm_dem_sniff_size(struct bmm_dem const* const dem,
    enum bmm_msg_num const num) {
  switch (num) {
    case BMM_MSG_NUM_ISTEP:
      return sizeof dem->time;
    case BMM_MSG_NUM_OPTS:
      return sizeof dem->opts;
    case BMM_MSG_NUM_NEIGH:
      return sizeof dem->cache + sizeof dem->link;
    case BMM_MSG_NUM_PARTS:
      return sizeof dem->part.n + sizeof dem->part;
  }

  dynamic_assert(false, "Unsupported message number");
}

bool bmm_dem_puts_stuff(struct bmm_dem const* const dem,
    enum bmm_msg_num const num) {
  switch (num) {
    case BMM_MSG_NUM_ISTEP:
      return msg_write(&dem->time, sizeof dem->time, NULL);
    case BMM_MSG_NUM_OPTS:
      return msg_write(&dem->opts, sizeof dem->opts, NULL);
    case BMM_MSG_NUM_NEIGH:
      return msg_write(&dem->cache, sizeof dem->cache, NULL) &&
        msg_write(&dem->link, sizeof dem->link, NULL);
    case BMM_MSG_NUM_PARTS:
      return msg_write(&dem->part.n, sizeof dem->part.n, NULL) &&
        msg_write(&dem->part, sizeof dem->part, NULL);
  }

  dynamic_assert(false, "Unsupported message number");
}

bool bmm_dem_puts(struct bmm_dem const* const dem,
    enum bmm_msg_num const num) {
  struct bmm_msg_spec spec;
  bmm_msg_spec_def(&spec);
  spec.msg.size = bmm_dem_sniff_size(dem, num) + BMM_MSG_NUMSIZE;

  return bmm_msg_spec_write(&spec, msg_write, NULL) &&
    bmm_msg_num_write(&num, msg_write, NULL) &&
    bmm_dem_puts_stuff(dem, num);
}

bool bmm_dem_cache_expired(struct bmm_dem const* const dem) {
  // TODO Use `dem->opts.part.rnew[1]` instead of `dem->part.r[ipart]`.
  double const r = dem->opts.cache.rcutoff / 2.0;

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    if (bmm_geom2d_cpdist2(dem->part.x[ipart], dem->cache.x[ipart],
          dem->opts.box.x, dem->opts.box.per) >=
        bmm_fp_sq(r - dem->part.r[ipart]))
      return true;

  return false;
}

__attribute__ ((__nonnull__))
static void bmm_dem_script_balance(struct bmm_dem* const dem) {
  double xcenter[BMM_NDIM];
  bmm_dem_est_center(xcenter, dem);

  double xcom[BMM_NDIM];
  bmm_dem_est_com(xcom, dem);

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    for (size_t idim = 0; idim < BMM_NDIM; ++idim) {
      dem->part.x[ipart][idim] += xcenter[idim] - xcom[idim];

      if (dem->opts.box.per[idim])
        dem->part.x[ipart][idim] = bmm_fp_uwrap(dem->part.x[ipart][idim],
            dem->opts.box.x[idim]);
    }
}

__attribute__ ((__nonnull__))
static bool bmm_dem_script_create_hc(struct bmm_dem* const dem) {
  double const etahc = bmm_geom_ballvol(0.5, BMM_NDIM);
  double const vhc = bmm_fp_prod(dem->opts.box.x, BMM_NDIM);
  double const eta = dem->opts.script.params[dem->script.i].create.eta;
  double const v = vhc * (etahc / eta);

  double const vlim = vhc * eta;
  // double const rspace = dem->opts.part.rnew[0];
  double const rspace = dem->opts.part.rnew[1];

  double x[BMM_NDIM];
  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    x[idim] = 0.0;

  double vnow = 0.0;

  for ever {
    double const r = bmm_random_get(dem->rng, dem->opts.part.rnew);
    double const v = bmm_geom_ballvol(r, BMM_NDIM);

    double const vnext = vnow + v;

    if (vnext < vlim) {
      size_t ipart = bmm_dem_addpart(dem);
      dem->part.r[ipart] = r;
      dem->part.m[ipart] = dem->opts.part.rho * v;

      // TODO Warning: nontermination.
      if (x[0] + 2.0 * r >= dem->opts.box.x[0]) {
        x[0] = 0.0;
        x[1] += 2.0 * rspace;
      }

      for (size_t idim = 0; idim < BMM_NDIM; ++idim)
        dem->part.x[ipart][idim] = x[idim] + dem->opts.part.rnew[1];

      // TODO Wow, lewd.
      for (size_t idim = 0; idim < BMM_NDIM; ++idim)
        dem->part.x[ipart][idim] += bmm_random_get(dem->rng,
            dem->opts.part.rnew) / 4.0;

      x[0] += 2.0 * rspace;

      vnow = vnext;
    } else
      break;
  }

  return true;
}

__attribute__ ((__nonnull__))
static bool bmm_dem_script_create_ball(struct bmm_dem* const dem) {
  return false;
}

__attribute__ ((__nonnull__))
static bool bmm_dem_script_create_gas(struct bmm_dem* const dem) {
  for (size_t ipart = 0; ipart < 64; ++ipart) {
    size_t const jpart = bmm_dem_addpart(dem);

    dem->part.r[jpart] = 0.03;
    dem->part.m[jpart] = 1.0;

    for (size_t idim = 0; idim < BMM_NDIM; ++idim)
      dem->part.x[jpart][idim] += gsl_rng_uniform(dem->rng) * dem->opts.box.x[idim];

      dem->part.omega[jpart] += (double) (rand() % 512 - 256);
  }
}

__attribute__ ((__nonnull__))
static bool bmm_dem_script_create_couple(struct bmm_dem* const dem) {
  size_t jpart;
  jpart = bmm_dem_addpart(dem);
  dem->part.r[jpart] = 0.03;
  dem->part.m[jpart] = 1.0;
  dem->part.x[jpart][0] += 0.45;
  dem->part.x[jpart][1] += 0.25;
  dem->part.v[jpart][0] += 1.0;
  dem->part.omega[jpart] += 400.0;
  jpart = bmm_dem_addpart(dem);
  dem->part.r[jpart] = 0.03;
  dem->part.m[jpart] = 1.0;
  dem->part.x[jpart][0] += 0.55;
  dem->part.x[jpart][1] += 0.25;
  dem->part.v[jpart][0] -= 1.0;
  dem->part.omega[jpart] += 400.0;
}

/// The call `bmm_dem_step(dem)`
/// advances the simulation `dem` by one step.
/// Make sure the simulation has not ended prior to the call
/// by calling `bmm_dem_script_ongoing` or `bmm_dem_script_trans`.
bool bmm_dem_step(struct bmm_dem* const dem) {
  switch (dem->opts.script.mode[dem->script.i]) {
    case BMM_DEM_MODE_CREATE:
      if (!bmm_dem_script_create_hc(dem))
        return false;

      bmm_dem_script_balance(dem);

      break;
    case BMM_DEM_MODE_TEST_GAS:
      if (!bmm_dem_script_create_gas(dem))
        return false;

      bmm_dem_script_balance(dem);

      break;
    case BMM_DEM_MODE_TEST_COUPLE:
      if (!bmm_dem_script_create_couple(dem))
        return false;

      bmm_dem_script_balance(dem);

      break;
    case BMM_DEM_MODE_SEDIMENT:
      // TODO This scheme sucks.
      dem->fext = BMM_DEM_FEXT_ABS;
      dem->ext.params.abs.fcohes = 1.0e+1;

      dem->fext = BMM_DEM_FEXT_HARM;
      dem->ext.params.harm.kcohes =
        dem->opts.script.params[dem->script.i].sediment.kcohes;

      break;
    case BMM_DEM_MODE_LINK:
      dem->fext = BMM_DEM_FEXT_NONE;

      if (!bmm_dem_link(dem))
        return false;

      break;
    case BMM_DEM_MODE_FAULT:

      break;
    case BMM_DEM_MODE_SEPARATE:

      break;
    case BMM_DEM_MODE_CRUNCH:

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
  bmm_dem_force(dem);
  bmm_dem_correct(dem);

  if (dem->time.istep % dem->opts.time.istab == 0)
    bmm_dem_stab(dem);

  dem->time.t += dem->opts.script.dt[dem->script.i];
  ++dem->time.istep;

  return true;
}

// TODO This looks just like `bmm_dem_script_trans`.
bool bmm_dem_comm(struct bmm_dem* const dem) {
  // TODO Make a mechanism to automate retransmission of differences only.

  double const toff = dem->time.t - dem->comm.tprev - dem->opts.comm.dt;

  if (toff >= 0.0) {
    dem->comm.tprev = dem->time.t;

    if (!bmm_dem_puts(dem, BMM_MSG_NUM_ISTEP))
      return false;

    if (!bmm_dem_puts(dem, BMM_MSG_NUM_NEIGH))
      return false;

    if (!bmm_dem_puts(dem, BMM_MSG_NUM_PARTS))
      return false;
  }

  return true;
}

static double abserr(double const x, double const z,
    __attribute__ ((__unused__)) void* const ptr) {
  return fabs(x) + z;
}

bool bmm_dem_report(struct bmm_dem const* const dem) {
  if (dem->opts.verbose) {
    if (fprintf(stderr, "Time Error: %g\n",
          bmm_fp_lfold(abserr,
            dem->script.toff, dem->opts.script.n, 0.0, NULL)) < 0)
      return false;
  }

  return true;
}

static bool bmm_dem_run_(struct bmm_dem* const dem) {
  int const sigs[] = {SIGINT, SIGQUIT, SIGTERM, SIGPIPE};
  if (bmm_sig_register(sigs, nmembof(sigs)) != SIZE_MAX) {
    BMM_TLE_STDS();

    return false;
  }

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

bool bmm_dem_run(struct bmm_dem* const dem) {
  bool const run = bmm_dem_run_(dem);
  bool const report = bmm_dem_report(dem);

  return run && report;
}

static bool bmm_dem_run_with_(struct bmm_dem* const dem) {
  gsl_rng_type const* const t = gsl_rng_env_setup();
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

bool bmm_dem_run_with(struct bmm_dem_opts const* const opts) {
  struct bmm_dem* const dem = malloc(sizeof *dem);
  if (dem == NULL) {
    BMM_TLE_STDS();

    return false;
  }

  bmm_dem_def(dem, opts);

  bool const result = bmm_dem_run_with_(dem);

  free(dem);

  return result;
}
