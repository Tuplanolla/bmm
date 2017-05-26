#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
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
#include "moore.h"
#include "msg.h"
#include "sig.h"
#include "size.h"
#include "tle.h"

// TODO Use `dynamic_assert(n >= 3, "Too few neighbor cells")` as appropriate.

void bmm_dem_ijcell(size_t* const pijcell,
    struct bmm_dem const* const dem, size_t const ipart) {
  for (size_t idim = 0; idim < BMM_NDIM; ++idim) {
    size_t const n = dem->opts.cache.ncell[idim];

    size_t const j = n - 1;
    double const r = (double) j;

    double const x = bmm_fp_lerp(dem->part.x[ipart][idim],
        0.0, dem->opts.box.x[idim], 1.0, r);

    if (x < 1.0)
      pijcell[idim] = 0;
    else if (x >= r)
      pijcell[idim] = j;
    else {
      size_t const k = (size_t) x;

      dynamic_assert(k < j, "Invalid truncation");

      pijcell[idim] = k;
    }
  }
}

// TODO Think about this.
// Relying on index space inequalities makes the cache unstable,
// since the addition and removal of particles permutes it.
bool bmm_dem_isneigh(struct bmm_dem* const dem,
    size_t const ipart, size_t const jpart, bool const lex) {
  if (ipart == jpart)
    return false;

  if (lex && ipart > jpart)
    return false;

  if (bmm_geom2d_cpdist2(dem->part.x[ipart], dem->part.x[jpart],
        dem->opts.box.x, dem->opts.box.per) >
      bmm_fp_sq(dem->opts.cache.rcutoff))
    return false;

  return true;
}

// TODO Name these and share them with the refresh mechanism.

// Usually you first fill `cache.part` by going over all particles.
// Then you go over all particles again and,
// for each particle in the same cell or upper Moore half,
// check the neighborhood condition.
//
// In this case the beginning is the same.
// However you only go over this particle and,
// for each particle in the same cell or upper Moore half,
// check the neighborhood condition.
// Additionally you also go through every particle in the reduced lower half,
// checking the neighborhood condition only for this particle.
//
// Particle removal is also a bit tricky since the last particle
// may suddenly change its index to some smaller value,
// which invalidates the cache of every particle in its cell wrt it.

/// The call `bmm_dem_cache_cell(dem, ipart)`
/// adds the new particle `ipart` to the appropriate neighbor cell
/// unless the cell group is full.
bool bmm_dem_cache_cell(struct bmm_dem* const dem, size_t const ipart) {
  bmm_dem_ijcell(dem->cache.cell[ipart], dem, ipart);

  size_t const icell = bmm_size_unhcd(dem->cache.cell[ipart],
      BMM_NDIM, dem->opts.cache.ncell);

  if (dem->cache.part[icell].n >= nmembof(dem->cache.part[icell].i))
    return false;

  dem->cache.part[icell].i[dem->cache.part[icell].n] = ipart;

  ++dem->cache.part[icell].n;

  return true;
}

/// The call `bmm_dem_cache_recell(dem)`
/// moves every particle to the appropriate neighbor cell
/// except for those particles whose cell group becomes full.
bool bmm_dem_cache_recell(struct bmm_dem* const dem) {
  for (size_t icell = 0; icell < nmembof(dem->cache.part); ++icell)
    dem->cache.part[icell].n = 0;

  bool p = true;

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    if (!bmm_dem_cache_cell(dem, ipart))
      p = false;

  return p;
}

// TODO Make this reduced and split off the special case of the center.

// Goes over these neighbor cells and adds them to `ipart`'s neighbors.
//
//       ^
//       | - + +
//     y | - + +
//       | - - +
//       +------->
//           x
//
bool bmm_dem_cache_doto(struct bmm_dem* const dem, size_t const ipart) {
  size_t ij[BMM_NDIM];

  size_t const nmoore = bmm_moore_ncpuh(ij,
      dem->cache.cell[ipart], BMM_NDIM,
      dem->opts.cache.ncell, dem->opts.box.per);

  for (size_t imoore = 0; imoore < nmoore; ++imoore) {
    bmm_moore_ijcpuh(ij, dem->cache.cell[ipart], imoore,
        BMM_NDIM, dem->opts.cache.ncell, dem->opts.box.per);

    size_t const icell = bmm_size_unhcd(ij, BMM_NDIM, dem->opts.cache.ncell);

    for (size_t igroup = 0; igroup < dem->cache.part[icell].n; ++igroup) {
      size_t const jpart = dem->cache.part[icell].i[igroup];

      if (bmm_dem_isneigh(dem, ipart, jpart, imoore == 0)) {
        if (dem->cache.neigh[ipart].n >= nmembof(dem->cache.neigh[ipart].i))
          return false;

        dem->cache.neigh[ipart].i[dem->cache.neigh[ipart].n] = jpart;

        ++dem->cache.neigh[ipart].n;
      }
    }
  }

  return true;
}

// Goes over these neighbor cells and adds `ipart` to their neighbors.
//
//       ^
//       | + - -
//     y | + - -
//       | + + -
//       +------->
//           x
//
bool bmm_dem_cache_dofrom(struct bmm_dem* const dem, size_t const ipart) {
  bool p = true;

  size_t ij[BMM_NDIM];

  size_t const nmoore = bmm_moore_ncplhr(ij,
      dem->cache.cell[ipart], BMM_NDIM,
      dem->opts.cache.ncell, dem->opts.box.per);

  for (size_t imoore = 0; imoore < nmoore; ++imoore) {
    bmm_moore_ijcplhr(ij, dem->cache.cell[ipart], imoore,
        BMM_NDIM, dem->opts.cache.ncell, dem->opts.box.per);

    size_t const icell = bmm_size_unhcd(ij, BMM_NDIM, dem->opts.cache.ncell);

    for (size_t igroup = 0; igroup < dem->cache.part[icell].n; ++igroup) {
      size_t const jpart = dem->cache.part[icell].i[igroup];

      if (bmm_dem_isneigh(dem, ipart, jpart, false)) {
        if (dem->cache.neigh[jpart].n >= nmembof(dem->cache.neigh[jpart].i))
          p = false;

        dem->cache.neigh[jpart].i[dem->cache.neigh[jpart].n] = ipart;

        ++dem->cache.neigh[jpart].n;
      }
    }
  }

  return p;
}

bool bmm_dem_reneigh(struct bmm_dem* const dem) {
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    dem->cache.neigh[ipart].n = 0;

    (void) bmm_dem_cache_doto(dem, ipart);
  }

  return true;
}

size_t bmm_dem_inspart(struct bmm_dem* const dem,
    double const r, double const m) {
  size_t const ipart = dem->part.n;

  if (ipart >= BMM_KPART)
    return BMM_KPART;

  ++dem->part.n;

  dem->part.l[ipart] = dem->part.lnew;
  ++dem->part.lnew;

  dem->part.role[ipart] = BMM_DEM_ROLE_FREE;

  dem->part.r[ipart] = r;
  dem->part.m[ipart] = m;
  dem->part.j[ipart] = bmm_geom_ballmoi(m, r, BMM_NDIM);

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

  // "In this case..."
  (void) bmm_dem_cache_cell(dem, ipart);

  // "However..."
  dem->cache.neigh[ipart].n = 0;

  (void) bmm_dem_cache_doto(dem, ipart);

  // "Additionally..."
  (void) bmm_dem_cache_dofrom(dem, ipart);

  return ipart;
}

bool bmm_dem_delpart(struct bmm_dem* const dem, size_t const ipart) {
  dynamic_assert(ipart < dem->part.n, "Index out of bounds");

  size_t const jpart = dem->part.n - 1;

  {
    dem->part.l[ipart] = dem->part.l[jpart];

    dem->part.role[ipart] = dem->part.role[jpart];

    dem->part.r[ipart] = dem->part.r[jpart];
    dem->part.m[ipart] = dem->part.m[jpart];
    dem->part.j[ipart] = dem->part.j[jpart];

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
  }

  {
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

  {
    for (size_t idim = 0; idim < BMM_NDIM; ++idim)
      dem->cache.cell[ipart][idim] = dem->cache.cell[jpart][idim];

    // The cell the particle being removed is registered to.
    size_t const icell = bmm_size_unhcd(dem->cache.cell[jpart],
        BMM_NDIM, dem->opts.cache.ncell);

    for (size_t kpart = 0; kpart < dem->cache.part[icell].n; ++kpart)
      if (dem->cache.part[icell].i[kpart] == ipart) {
        size_t const lpart = dem->cache.part[icell].n - 1;

        dem->cache.part[icell].i[kpart] = dem->cache.part[icell].i[lpart];

        --dem->cache.part[icell].n;

        // TODO What if not found?
        break;
      }

    dem->cache.neigh[ipart].n = dem->cache.neigh[jpart].n;

    for (size_t kpart = 0; kpart < dem->cache.neigh[jpart].n; ++kpart)
      dem->cache.neigh[ipart].i[kpart] = dem->cache.neigh[jpart].i[kpart];

    dem->cache.neigh[ipart].n = dem->cache.neigh[jpart].n;

    for (size_t kpart = 0; kpart < dem->cache.neigh[jpart].n; ++kpart)
      dem->cache.neigh[ipart].i[kpart] = dem->cache.neigh[jpart].i[kpart];
  }

  --dem->part.n;

  // TODO Factor index structure adjustment.
  // This is the slow part.
  {
    for (size_t kpart = 0; kpart < dem->part.n; ++kpart)
      for (size_t lpart = 0; lpart < dem->link.part[kpart].n; ++lpart)
        if (dem->link.part[kpart].i[lpart] == jpart) {
          dem->link.part[kpart].i[lpart] = ipart;

          break;
        }

    size_t const ncell = bmm_size_prod(dem->opts.cache.ncell, BMM_NDIM);

    for (size_t icell = 0; icell < ncell; ++icell)
      for (size_t kpart = 0; kpart < dem->cache.part[icell].n; ++kpart)
        if (dem->cache.part[icell].i[kpart] == jpart) {
          dem->cache.part[icell].i[kpart] = ipart;

          break;
        }

    for (size_t kpart = 0; kpart < dem->part.n; ++kpart)
      for (size_t lpart = 0; lpart < dem->cache.neigh[kpart].n; ++lpart)
        if (dem->cache.neigh[kpart].i[lpart] == jpart) {
          dem->cache.neigh[kpart].i[lpart] = ipart;

          break;
        }
  }

  return true;
}

typedef void (* bmm_fpair)(struct bmm_dem*, size_t, size_t);

void bmm_dem_force_pair(struct bmm_dem* const dem,
    size_t const ipart, size_t const jpart) {
  double xdiff[BMM_NDIM];
  bmm_geom2d_cpdiff(xdiff, dem->part.x[ipart], dem->part.x[jpart],
      dem->opts.box.x, dem->opts.box.per);

  double const d2 = bmm_geom2d_norm2(xdiff);

  double const r = dem->part.r[ipart] + dem->part.r[jpart];
  double const r2 = bmm_fp_sq(r);

  if (d2 > r2)
    return;

  // TODO Maybe cache some of these.

  double const d = sqrt(d2);

  double xnorm[BMM_NDIM];
  bmm_geom2d_scale(xnorm, xdiff, 1.0 / d);

  double xtang[BMM_NDIM];
  bmm_geom2d_rperp(xtang, xnorm);

  double vdiff[BMM_NDIM];
  bmm_geom2d_diff(vdiff, dem->part.v[ipart], dem->part.v[jpart]);

  // TODO Apply 3.2 for real now.
  double const xi = r - d;
  double const dotxi = bmm_geom2d_dot(vdiff, xnorm);
  double const vtang = bmm_geom2d_dot(vdiff, xtang) +
    dem->part.r[ipart] * dem->part.omega[ipart] +
    dem->part.r[jpart] * dem->part.omega[jpart];

  double fn = 0.0;
  switch (dem->opts.fnorm) {
    case BMM_DEM_FNORM_DASHPOT:
      fn = fmax(0.0, dem->opts.part.y * xi +
          dem->opts.norm.params.dashpot.gamma * dotxi);

      break;
  }

  // TODO Investigate the sign.
  double fdiff[2];
  bmm_geom2d_scale(fdiff, xnorm, -fn);

  bmm_geom2d_addto(dem->part.f[ipart], fdiff);

  double fdiff2[2];
  bmm_geom2d_scale(fdiff2, fdiff, -1.0);

  bmm_geom2d_addto(dem->part.f[jpart], fdiff2);

  double ft = 0.0;
  switch (dem->opts.ftang) {
    case BMM_DEM_FTANG_HW:
      ft = -copysign(fmin(dem->opts.tang.params.hw.gamma * fabs(vtang),
            dem->opts.tang.params.hw.mu * fn), vtang);

      break;
  }

  bmm_geom2d_scale(fdiff, xtang, ft);

  bmm_geom2d_addto(dem->part.f[ipart], fdiff);

  dem->part.tau[ipart] += ft * dem->part.r[ipart];

  bmm_geom2d_scale(fdiff2, fdiff, -1.0);

  bmm_geom2d_addto(dem->part.f[jpart], fdiff2);

  dem->part.tau[jpart] += ft * dem->part.r[jpart];
}

void bmm_dem_force_ambient(struct bmm_dem* const dem, size_t const ipart) {
  // TODO Shimmy these switches outside the loops,
  // perhaps by passing in function pointers.
  switch (dem->opts.famb) {
    case BMM_DEM_FAMB_CREEPING:
      for (size_t idim = 0; idim < nmembof(dem->part.f[ipart]); ++idim)
        dem->part.f[ipart][idim] *= 1.0;

      break;
  }
}

void bmm_dem_force_link(struct bmm_dem* const dem,
    size_t const ipart, size_t const jpart, size_t const ilink) {
  // TODO This.
  return;
}

void bmm_dem_force(struct bmm_dem* const dem) {
  // TODO Point these to functions.
  bmm_fpair fnorm = NULL;
  bmm_fpair ftang = NULL;

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    for (size_t idim = 0; idim < nmembof(dem->part.f[ipart]); ++idim)
      dem->part.f[ipart][idim] = 0.0;

    dem->part.tau[ipart] = 0.0;
  }

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    bmm_dem_force_ambient(dem, ipart);

  switch (dem->opts.caching) {
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

  // TODO Calculate force feedback from the residuals.
}

void bmm_dem_integ_euler(struct bmm_dem* const dem) {
  double const dt = dem->opts.script.dt[dem->script.i];

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    for (size_t idim = 0; idim < 2; ++idim) {
      dem->part.x[ipart][idim] = dem->part.x[ipart][idim] + dem->part.v[ipart][idim] * dt;
      if (dem->opts.box.per[idim])
        dem->part.x[ipart][idim] = bmm_fp_uwrap(dem->part.x[ipart][idim],
            dem->opts.box.x[idim]);

      // TODO Ordering.
      dem->part.v[ipart][idim] = dem->part.v[ipart][idim] + dem->part.a[ipart][idim] * dt;

      dem->part.a[ipart][idim] = dem->part.f[ipart][idim] / dem->part.m[ipart];
    }

    dem->part.phi[ipart] = dem->part.phi[ipart] + dem->part.omega[ipart] * dt;

    dem->part.omega[ipart] = dem->part.omega[ipart] + dem->part.alpha[ipart] * dt;

    dem->part.alpha[ipart] = dem->part.tau[ipart] / dem->part.j[ipart];
  }
}

/// Clamp periodic domains, renormalize normal vectors, etc.
void bmm_dem_stab(struct bmm_dem* const dem) {
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    dem->part.phi[ipart] = bmm_fp_uwrap(dem->part.phi[ipart], M_2PI);
}

void bmm_dem_predict(struct bmm_dem* const dem) {
  switch (dem->opts.integ) {
    case BMM_DEM_INTEG_EULER:
      break;
    case BMM_DEM_INTEG_GEAR:
      // bmm_dem_integ_gearpred(dem);

      break;
  }
}

void bmm_dem_correct(struct bmm_dem* const dem) {
  switch (dem->opts.integ) {
    case BMM_DEM_INTEG_EULER:
      bmm_dem_integ_euler(dem);

      break;
    case BMM_DEM_INTEG_GEAR:
      // bmm_dem_integ_gearcorr(dem);

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
  switch (dem->opts.flink) {
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
  switch (dem->opts.caching) {
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

// TODO These are dubious for empty sets.

// Maximum velocity estimator.
void bmm_dem_maxvel(double* const v, struct bmm_dem const* const dem) {
  for (size_t idim = 0; idim < 2; ++idim) {
    v[idim] = 0.0;

    for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
      v[idim] = fmax(v[idim], dem->part.v[ipart][idim]);
  }
}

// Maximum radius estimator.
double bmm_dem_maxrad(struct bmm_dem const* const dem) {
  double r = 0.0;

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    r = fmax(r, dem->part.r[ipart]);

  return r;
}

// Drift time estimator.
double bmm_dem_drift(struct bmm_dem const* const dem) {
  double t = (double) INFINITY;

  // TODO Make this a configuration constant.
  double const rad = bmm_dem_maxrad(dem);

  double v[2];
  bmm_dem_maxvel(v, dem);

  for (size_t idim = 0; idim < 2; ++idim)
    t = fmin(t, (0.5 *
          dem->opts.box.x[idim] / (double) dem->opts.cache.ncell[idim]) - rad) /
        (v[idim] + 0.01);

  return t;
}

// Total kinetic energy estimator.
double bmm_dem_ekinetic(struct bmm_dem const* const dem) {
  double e = 0.0;

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    for (size_t idim = 0; idim < 2; ++idim)
      e += dem->part.m[ipart] * bmm_fp_sq(dem->part.v[ipart][idim]);

    e += dem->part.j[ipart] * bmm_fp_sq(dem->part.omega[ipart]);
  }

  return e * 0.5;
}

// Total momentum estimator.
double bmm_dem_pvector(struct bmm_dem const* const dem) {
  double p[2];
  for (size_t idim = 0; idim < 2; ++idim)
    p[idim] = 0.0;

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    for (size_t idim = 0; idim < 2; ++idim)
      p[idim] += dem->part.m[ipart] * dem->part.v[ipart][idim];

  return bmm_geom2d_norm(p);
}

// Individual momentum estimator.
double bmm_dem_pscalar(struct bmm_dem const* const dem) {
  double p = 0.0;

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    p += dem->part.m[ipart] * bmm_geom2d_norm(dem->part.v[ipart]);

  return p;
}

// Mean coefficient of restitution
// (just linear dashpot for now, also a bit wrong).
double bmm_dem_cor(struct bmm_dem const* const dem) {
  double e = 0.0;

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    double const mred =
      dem->part.m[ipart] * dem->part.m[ipart] /
      (dem->part.m[ipart] + dem->part.m[ipart]);
    e += exp(-M_PI * dem->opts.norm.params.dashpot.gamma / (2.0 * mred) /
        sqrt(dem->opts.part.y / mred -
          bmm_fp_sq(dem->opts.norm.params.dashpot.gamma / (2.0 * mred))));
  }

  return e / (double) dem->part.n;
}

void bmm_dem_opts_def(struct bmm_dem_opts* const opts) {
  // This is here just to help Valgrind and cover up my mistakes.
  (void) memset(opts, 0, sizeof *opts);

  opts->verbose = true;

  opts->init = BMM_DEM_INIT_TRIAL;
  opts->integ = BMM_DEM_INTEG_EULER;
  opts->caching = BMM_DEM_CACHING_NONE;
  opts->caching = BMM_DEM_CACHING_NEIGH;
  opts->famb = BMM_DEM_FAMB_CREEPING;
  opts->fnorm = BMM_DEM_FNORM_DASHPOT;
  opts->ftang = BMM_DEM_FTANG_HW;

  opts->norm.params.dashpot.gamma = 1.0;
  opts->tang.params.hw.gamma = 1.0;
  opts->tang.params.hw.mu = 1.0;

  for (size_t idim = 0; idim < nmembof(opts->box.x); ++idim)
    opts->box.x[idim] = 1.0;

  for (size_t idim = 0; idim < nmembof(opts->box.per); ++idim)
    opts->box.per[idim] = false;

  opts->ambient.eta = 0.0;

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

  for (size_t idim = 0; idim < nmembof(opts->cache.ncell); ++idim)
    opts->cache.ncell[idim] = 5;

  opts->cache.rcutoff = 1.0;
}

void bmm_dem_def(struct bmm_dem* const dem,
    struct bmm_dem_opts const* const opts) {
  // This is here just to help Valgrind and cover up my mistakes.
  (void) memset(dem, 0, sizeof *dem);

  dem->opts = *opts;

  dem->time.istep = 0;
  dem->time.istab = 1000;
  dem->time.t = 0.0;

  dem->part.n = 0;
  dem->part.lnew = 0;

  for (size_t ipart = 0; ipart < nmembof(dem->link.part); ++ipart)
    dem->link.part[ipart].n = 0;

  dem->script.i = 0;
  dem->script.tprev = 0.0;

  dem->comm.lnew = 0;
  dem->comm.tprev = (double) NAN;
  dem->comm.tnext = 0.0;

  dem->cache.i = 0;
  dem->cache.tpart = (double) NAN;
  dem->cache.tprev = (double) NAN;
  dem->cache.tnext = 0.0;

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

bool bmm_dem_script_pushidle(struct bmm_dem_opts* const opts,
    double const tspan) {
  if (opts->script.n >= nmembof(opts->script.mode))
    return false;

  double const dt = 1e-3;
  enum bmm_dem_mode mode = BMM_DEM_MODE_IDLE;

  opts->script.tspan[opts->script.n] = tspan;
  opts->script.dt[opts->script.n] = dt;
  opts->script.mode[opts->script.n] = mode;

  ++opts->script.n;

  return true;
}

bool bmm_dem_script_done(struct bmm_dem const* const dem) {
  return dem->script.i >= dem->opts.script.n;
}

enum bmm_dem_step bmm_dem_script_trans(struct bmm_dem* const dem) {
  if (bmm_dem_script_done(dem))
    return BMM_DEM_STEP_STOP;

  double const toff = dem->time.t -
    (dem->script.tprev + dem->opts.script.tspan[dem->script.i]);

  if (toff >= 0.0) {
    dem->script.tprev = dem->time.t;
    dem->script.ttrans[dem->script.i] = dem->time.t;
    dem->script.toff[dem->script.i] = toff;

    ++dem->script.i;

    if (bmm_dem_script_done(dem))
      return BMM_DEM_STEP_STOP;
  }

  return BMM_DEM_STEP_CONT;
}

enum bmm_dem_step bmm_dem_step(struct bmm_dem* const dem) {
  switch (bmm_dem_script_trans(dem)) {
    case BMM_DEM_STEP_ERROR:
      return BMM_DEM_STEP_ERROR;
    case BMM_DEM_STEP_STOP:
      return BMM_DEM_STEP_STOP;
  }

  size_t const istage = dem->script.i;

  switch (dem->opts.script.mode[istage]) {
    case BMM_DEM_MODE_BEGIN:

      break;
    case BMM_DEM_MODE_SEDIMENT:

      break;
    case BMM_DEM_MODE_LINK:

      break;
    case BMM_DEM_MODE_SMASH:

      break;
    case BMM_DEM_MODE_ACCEL:

      break;
  }

  if (dem->time.t >= dem->cache.tnext) {
    bmm_dem_cache_recell(dem);

    bmm_dem_reneigh(dem);

    // TODO Take `tstep` into account
    // since the update only happens after the drift time has passed.
    // dem->cache.tnext += bmm_dem_drift(dem);
  }

  bmm_dem_predict(dem);

  bmm_dem_force(dem);

  bmm_dem_correct(dem);

  if (dem->time.istep % dem->time.istab == 0)
    bmm_dem_stab(dem);

  ++dem->time.istep;

  dem->time.t += dem->opts.script.dt[istage];

  return BMM_DEM_STEP_CONT;
}

// TODO This ought to be const-qualified.
bool bmm_dem_comm(struct bmm_dem* const dem) {
  // TODO Make a mechanism to automate retransmission of differences only.

  bmm_dem_puts(dem, BMM_MSG_NUM_ISTEP);
  bmm_dem_puts(dem, BMM_MSG_NUM_NEIGH);
  bmm_dem_puts(dem, BMM_MSG_NUM_PARTS);

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

// TODO Would it be beneficial to be higher-order?
static bool bmm_dem_run_(struct bmm_dem* const dem) {
  int const sigs[] = {SIGINT, SIGQUIT, SIGTERM, SIGPIPE};
  if (bmm_sig_register(sigs, nmembof(sigs)) != SIZE_MAX) {
    BMM_TLE_STDS();

    return false;
  }

  // TODO Get rid of these after refactoring `dem.c`.

  dem->opts.box.per[0] = true;
  dem->opts.box.per[1] = false;

  dem->opts.cache.rcutoff = 1.0;

  dem->opts.part.y = 1e+4;

  bmm_dem_script_pushidle(&dem->opts, 0.04);
  bmm_dem_script_pushidle(&dem->opts, 0.26);

  // Random stuff.
  for (size_t ipart = 0; ipart < 64; ++ipart) {
    size_t const jpart = bmm_dem_inspart(dem, 0.03, 1.0);

    for (size_t idim = 0; idim < nmembof(dem->part.x[jpart]); ++idim)
      dem->part.x[jpart][idim] += gsl_rng_uniform(dem->rng) * dem->opts.box.x[idim];

      dem->part.omega[jpart] += (double) (rand() % 512 - 256);
  }

  // Rotating couple.
  size_t jpart;
  jpart = bmm_dem_inspart(dem, 0.03, 1.0);
  dem->part.x[jpart][0] += 0.45;
  dem->part.x[jpart][1] += 0.5;
  dem->part.v[jpart][0] += 1.0;
  dem->part.omega[jpart] += 400.0;
  jpart = bmm_dem_inspart(dem, 0.03, 1.0);
  dem->part.x[jpart][0] += 0.55;
  dem->part.x[jpart][1] += 0.5;
  dem->part.v[jpart][0] -= 1.0;
  dem->part.omega[jpart] += 400.0;

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

    if (!bmm_dem_comm(dem))
      return false;

    switch (bmm_dem_step(dem)) {
      case BMM_DEM_STEP_ERROR:
        return false;
      case BMM_DEM_STEP_STOP:
        return true;
    }
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
