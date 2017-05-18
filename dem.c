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

#ifdef _GNU_SOURCE
#ifdef DEBUG
#include <fenv.h>
#endif
#endif

// TODO Wow, disgusting.

/*
extern inline void bmm_dem_clear(struct bmm_dem_list* const list);

extern inline bool bmm_dem_push(struct bmm_dem_list* const list, size_t const x);

extern inline size_t bmm_dem_size(struct bmm_dem_list const* const list);

extern inline size_t bmm_dem_get(struct bmm_dem_list const* const list, size_t const i);
*/

void bmm_dem_ijcell(size_t* const pijcell,
    struct bmm_dem const* const dem, size_t const ipart) {
  for (size_t idim = 0; idim < BMM_NDIM; ++idim) {
    size_t const n = dem->opts.cache.ncell[idim];

    // TODO How about periodicity?
    dynamic_assert(n >= 3, "Too few neighbor cells");

    size_t const j = n - 1;
    double const r = (double) j - 1.0;
    double const x = bmm_fp_lerp(dem->part.x[ipart][idim],
        0.0, dem->opts.box.x[idim], 0.0, r);

    if (x < 0.0)
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

bool bmm_dem_isneigh(struct bmm_dem* const dem,
    size_t const ipart, size_t const jpart) {
  return bmm_geom2d_cpdist2(dem->part.x[ipart], dem->part.x[jpart],
      dem->opts.box.x, dem->opts.box.per) <=
    bmm_fp_sq(dem->opts.cache.rcutoff);
}

// TODO Consider a sensible implementation.
#define bmm_dem_push(p, x) \
  begin \
    if (p.n >= nmembof(p.i)) \
      goto br; \
    p.i[p.n] = x; \
    ++p.n; \
  end

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

/// The call `bmm_dem_cache_cell(dem, ipart)`
/// adds the new particle `ipart` to the appropriate neighbor cell
/// unless the cell group is full.
bool bmm_dem_cache_cell(struct bmm_dem* const dem, size_t const ipart) {
  bmm_dem_ijcell(dem->cache.cell[ipart], dem, ipart);

  size_t const icell = bmm_size_unhc(dem->cache.cell[ipart],
      BMM_NDIM, BMM_NCELL);

  bmm_dem_push(dem->cache.part[icell], ipart);

  return true;

br:
  return false;
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

// Goes over these neighbor cells and adds them to `ipart`'s neighbors.
//
//       ^
//       | - + +
//     y | - + +
//       | - - +
//       +------->
//           x
//
bool bmm_dem_for2(struct bmm_dem* const dem, size_t const ipart) {
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

      if (bmm_dem_isneigh(dem, ipart, jpart))
        bmm_dem_push(dem->cache.neigh[ipart], jpart);
    }
  }

  return true;

br:
  return false;
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
bool bmm_dem_for3(struct bmm_dem* const dem, size_t const ipart) {
  bool p = true;

  size_t ij[BMM_NDIM];

  size_t const nmoore = bmm_moore_ncplhr(ij,
      dem->cache.cell[ipart], BMM_NDIM,
      dem->opts.cache.ncell, dem->opts.box.per);

  for (size_t imoore = 0; imoore < nmoore; ++imoore) {
    bmm_moore_ijcpuh(ij, dem->cache.cell[ipart], imoore,
        BMM_NDIM, dem->opts.cache.ncell, dem->opts.box.per);

    size_t const icell = bmm_size_unhcd(ij, BMM_NDIM, dem->opts.cache.ncell);

    for (size_t igroup = 0; igroup < dem->cache.part[icell].n; ++igroup) {
      size_t const jpart = dem->cache.part[icell].i[igroup];

      if (bmm_dem_isneigh(dem, ipart, jpart))
        bmm_dem_push(dem->cache.neigh[jpart], ipart);
    }

    continue;

br:
    p = false;
  }

  return p;
}

size_t bmm_dem_addpart(struct bmm_dem* const dem,
    double const r, double const m) {
  size_t const ipart = dem->part.n;

  if (ipart >= BMM_NPART)
    return BMM_NPART;

  ++dem->part.n;

  dem->part.l[ipart] = dem->part.lnew;
  ++dem->part.lnew;

  dem->part.role[ipart] = BMM_DEM_ROLE_FREE;

  dem->part.r[ipart] = r;
  dem->part.m[ipart] = m;
  dem->part.j[ipart] = m * bmm_geom_ballmoi(r, BMM_NDIM);

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

  (void) bmm_dem_for2(dem, ipart);

  // "Additionally..."
  (void) bmm_dem_for3(dem, ipart);

  return ipart;
}

bool bmm_dem_rmpart(struct bmm_dem* const dem, size_t const ipart) {
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
      for (size_t iend = 0; iend < 2; ++iend)
        dem->link.part[ipart].phirest[ilink][iend] =
          dem->link.part[jpart].phirest[ilink][iend];

    for (size_t ilink = 0; ilink < dem->link.part[jpart].n; ++ilink)
      dem->link.part[ipart].ftens[ilink] = dem->link.part[jpart].ftens[ilink];

    for (size_t ilink = 0; ilink < dem->link.part[jpart].n; ++ilink)
      dem->link.part[ipart].fshear[ilink] = dem->link.part[jpart].fshear[ilink];
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
  // TODO Make a predicate function.
  double dr[2];
  bmm_geom2d_pdiff(dr,
      dem->part.x[ipart],
      dem->part.x[jpart],
      dem->opts.box.x);

  double const d2 = bmm_geom2d_norm2(dr);

  double const r = dem->part.r[ipart] + dem->part.r[jpart];

  double const r2 = bmm_fp_sq(r);

  if (d2 < r2) {
    double const x = r - sqrt(d2);

    // TODO Write a fast vector shift function.
    listy->thingy[ineigh].x[1] = listy->thingy[ineigh].x[0];
    listy->thingy[ineigh].x[0] = x;

    // TODO Write a thing for this too.
    double const dx = listy->thingy[ineigh].x[0] - listy->thingy[ineigh].x[1];
    double const v = dx / dem->opts.tstep;

    // TODO Use getters.
    double const f = fmax(0.0,
        dem->opts.part.y * x + dem->opts.yelast * v);

    double n[2];
    bmm_geom2d_normal(n, dr);

    double t[2];
    bmm_geom2d_rperp(t, n);

    double df[2];
    bmm_geom2d_scale(df, n, -f);

    bmm_geom2d_add(dem->part.f[ipart], dem->part.f[ipart], df);
  }
}

void bmm_dem_force_link(struct bmm_dem* const dem,
    size_t const ipart, size_t const jpart, size_t const ilink) {
  // TODO Copy-pasted from above...

  double dr[2];
  bmm_geom2d_pdiff(dr,
      dem->part.x[ipart],
      dem->part.x[jpart],
      dem->opts.box.x);

  double const d = bmm_geom2d_norm(dr);

  // TODO Use getters.
  double const x0 = list->linkl[ilink].x0;
  double const f = dem->opts.klink * (x0 - d);

  double n[2];
  bmm_geom2d_normal(n, dr);

  double t[2];
  bmm_geom2d_rperp(t, n);

  double df[2];
  bmm_geom2d_scale(df, n, -f);

  bmm_geom2d_add(dem->part.f[ipart], dem->part.f[ipart], df);
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

void bmm_dem_force(struct bmm_dem* const dem) {
  // TODO Point these to functions.
  bmm_fpair fnorm = NULL;
  bmm_fpair ftang = NULL;

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    for (size_t idim = 0; idim < nmembof(dem->part.f[ipart]); ++idim)
      dem->part.f[ipart][idim] = 0.0;

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
  // TODO Script system!
  double const dt = dem->opts.tstep;

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    for (size_t idim = 0; idim < 2; ++idim) {
      dem->part.a[ipart][idim] = dem->part.f[ipart][idim] / dem->part.m[ipart];

      dem->part.x[ipart][idim] = bmm_fp_uwrap(
          dem->part.x[ipart][idim] +
          dem->part.v[ipart][idim] * dt, dem->opts.box.x[idim]);

      dem->part.v[ipart][idim] = dem->part.v[ipart][idim] + dem->part.a[ipart][idim] * dt;
    }

    dem->part.alpha[ipart] = dem->part.tau[ipart] / dem->part.moi[ipart];

    dem->part.phi[ipart] = dem->part.phi[ipart] + dem->part.omega[ipart] * dt;

    dem->part.omega[ipart] = dem->part.omega[ipart] + dem->part.alpha[ipart] * dt;
  }
}

/// Clamp domains, renormalize normal vectors, etc.
void bmm_dem_stab(struct bmm_dem* const dem) {
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    dem->part.phi[ipart] = bmm_fp_uwrap(dem->part.phi[ipart], M_2PI);
}

void bmm_dem_predict(struct bmm_dem* const dem) {
  switch (dem->integ) {
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

// TODO Real physics.
bool bmm_dem_relink(struct bmm_dem* const dem) {
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    struct bmm_dem_listl* const list = &dem->buf.links[ipart];

    for (size_t ilink = 0; ilink < bmm_dem_sizel(list); ++ilink) {
      size_t const jpart = bmm_dem_getl(list, ilink);

      // TODO Copy-pasted from elsewhere...

      double dr[2];
      bmm_geom2d_pdiff(dr,
          dem->part.x[ipart],
          dem->part.x[jpart],
          dem->opts.box.x);

      double const d = bmm_geom2d_norm(dr);

      // TODO Use getters.
      double const x0 = list->linkl[ilink].x0;

      if (d > x0 * 1.75) {
        // TODO Remove element properly.
        --list->n;
        list->linkl[ilink] = list->linkl[list->n];
        --ilink;
      }
    }
  }

  return true;
}

// TODO Check triangulation quality and compare with Delaunay.
bool bmm_dem_link(struct bmm_dem* const dem) {
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    struct bmm_dem_listy* const listy = &dem->buf.neigh.neighs[ipart];

    struct bmm_dem_listl* const list = &dem->buf.links[ipart];

    bmm_dem_clearl(list);

    for (size_t ineigh = 0; ineigh < bmm_dem_sizey(listy); ++ineigh) {
      size_t const jpart = bmm_dem_gety(listy, ineigh);

      if (jpart == ipart)
        continue;

      double const d2 = bmm_geom2d_pdist2(
          dem->part.x[ipart],
          dem->part.x[jpart],
          dem->opts.box.x);

      if (d2 < bmm_fp_sq(dem->opts.linkslurp *
            (dem->part.r[ipart] + dem->part.r[jpart])))
        if (bmm_dem_pushl(list, jpart))
          list->linkl[list->n - 1].x0 = dem->opts.linkoff * sqrt(d2);
    }
  }

  return true;
}

// TODO All of it.
bool bmm_dem_break(struct bmm_dem* const dem) {
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    struct bmm_dem_listl* const list = &dem->buf.links[ipart];

    for (size_t ilink = 0; ilink < bmm_dem_sizel(list); ++ilink) {
      size_t const jpart = bmm_dem_getl(list, ilink);

      if (fabs(dem->part.x[ipart][1] - dem->opts.box.x[1] * 0.5) <
          0.1 * dem->opts.box.x[1]) {
        double p = gsl_rng_uniform(dem->rng);

        if (p > 0.1) {
          // TODO Copy-pasted from somewhere...
          --list->n;
          list->linkl[ilink] = list->linkl[list->n];
          --ilink;
        }
      }
    }

    // TODO This requires fundamental rebuilding of indices,
    // which deserves to be done separately.
    // This also causes nonphysical symmetry breaking.
    if (fabs(dem->part.x[ipart][1] - dem->opts.box.x[1] * 0.5) <
        0.1 * dem->opts.box.x[1]) {
      double p = gsl_rng_uniform(dem->rng);

      if (p > 1.0) {
        // TODO Copy-pasted from up there...
        --dem->part.n;
        dem->buf.parts[ipart] = dem->buf.parts[dem->part.n];
        dem->buf.partcs[ipart] = dem->buf.partcs[dem->part.n];
        dem->buf.neigh.neighs[ipart] = dem->buf.neigh.neighs[dem->part.n];
        dem->buf.links[ipart] = dem->buf.links[dem->part.n];
        --ipart;
      }
    }
  }

  return true;
}

bool bmm_dem_fix(struct bmm_dem* const dem) {
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    if (dem->part.x[ipart][1] > 0.8 * dem->opts.box.x[1]) {
      dem->part.free[ipart] = false;
      dem->part.nondr[ipart] = false;
    } else if (dem->part.x[ipart][1] < 0.2 * dem->opts.box.x[1])
      dem->part.free[ipart] = false;
  }

  return true;
}

bool bmm_dem_fault(struct bmm_dem* const dem) {
  for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
    if (dem->part.x[ipart][1] > 0.5 * dem->opts.box.x[1])
      dem->part.x[ipart][1] += dem->opts.yoink * dem->opts.box.x[1];
    else
      dem->part.x[ipart][1] -= dem->opts.yoink * dem->opts.box.x[1];
  }

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
    t = fmin(t, (fmin(
            dem->opts.rmax,
            0.5 * dem->opts.box.x[idim] / (double) dem->opts.ncell[idim]) - rad) /
        (v[idim] + dem->opts.vleeway));

  return t;
}

// Total kinetic energy estimator.
double bmm_dem_ekinetic(struct bmm_dem const* const dem) {
  double e = 0.0;

  for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
    for (size_t idim = 0; idim < 2; ++idim)
      e += dem->part.m[ipart] * bmm_fp_sq(dem->part.v[ipart][idim]);

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
    e += exp(-M_PI * dem->opts.yelast / (2.0 * mred) /
        sqrt(dem->opts.part.y / mred -
          bmm_fp_sq(dem->opts.yelast / (2.0 * mred))));
  }

  return e / (double) dem->part.n;
}

void bmm_dem_opts_def(struct bmm_dem_opts* const opts) {
  // This is here just to help Valgrind and cover up my mistakes.
  (void) memset(opts, 0, sizeof *opts);

  opts->init = BMM_DEM_INIT_TRIAL;
  opts->integ = BMM_DEM_INTEG_EULER;
  opts->caching = BMM_DEM_CACHING_NEIGH;
  opts->famb = BMM_DEM_FAMB_CREEPING;
  opts->fnorm = BMM_DEM_FNORM_DASHPOT;
  opts->ftang = BMM_DEM_FTANG_NONE;

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

  opts->script.n = 0;

  opts->comm.dt = 1.0;
  opts->comm.flip = true;
  opts->comm.flop = true;
  opts->comm.flap = true;
  opts->comm.flup = true;

  for (size_t idim = 0; idim < nmembof(opts->cache.ncell); ++idim)
    opts->cache.ncell[idim] = 3;

  opts->cache.rcutoff = 1.0;
}

// TODO Rewrite these.

static void bmm_bottom(struct bmm_dem* const dem) {
  double x = 0.0;

  for ever {
    double r = dem->opts.rmean + gsl_ran_gaussian(dem->rng, dem->opts.rsd);

    x += r;

    struct bmm_dem_partc* const partc = &dem->buf.partcs[dem->part.n];
    struct bmm_dem_part* const part = &dem->buf.parts[dem->part.n];

    bmm_dem_partc_def(partc);
    bmm_dem_part_def(part);

    partc->r = r;
    partc->free = false;
    part->lin.x[0] = x;
    part->lin.x[1] = dem->opts.box.x[1] / 32.0;

    if (x + r >= dem->opts.box.x[0]) {
      x -= r;

      double const rprime = (dem->opts.box.x[0] - x) / 2.0;

      x += rprime;

      partc->r = rprime;
      partc->free = false;
      part->lin.x[0] = x;
      part->lin.x[1] = dem->opts.box.x[1] / 32.0;

      ++dem->part.n;

      break;
    }

    x += r;

    ++dem->part.n;
  }
}

static bool bmm_disperse(struct bmm_dem* const dem) {
  size_t success = 0;

  size_t maxfail = 100;
  size_t fail = 0;

  while (fail < maxfail && dem->part.n < BMM_NPART) {
    double const r = dem->opts.rmean +
      gsl_ran_gaussian(dem->rng, dem->opts.rsd);

    double x[2];
    x[0] = gsl_rng_uniform(dem->rng) * dem->opts.box.x[0];
    x[1] = gsl_rng_uniform(dem->rng) * dem->opts.box.x[1];

    for (size_t ipart = 0; ipart < dem->part.n; ++ipart) {
      double dr[2];
      bmm_geom2d_pdiff(dr,
          dem->part.x[ipart],
          x,
          dem->opts.box.x);

      double const d2 = bmm_geom2d_norm2(dr);

      double const rs = dem->part.r[ipart] + r;

      double const r2 = bmm_fp_sq(rs);

      if (d2 < r2) {
        ++fail;
        goto ret;
      }
    }

    fail = 0;

    struct bmm_dem_partc* const partc = &dem->buf.partcs[dem->part.n];
    struct bmm_dem_part* const part = &dem->buf.parts[dem->part.n];

    bmm_dem_partc_def(partc);
    bmm_dem_part_def(part);

    partc->r = r;
    partc->free = true;
    part->lin.x[0] = x[0];
    part->lin.x[1] = x[1];

    ++dem->part.n;

    ++success;
ret: ;
  }

  // TODO The lucky number.
  return success > dem->opts.lucky;
}

void bmm_dem_def(struct bmm_dem* const dem,
    struct bmm_dem_opts const* const opts) {
  // This is here just to help Valgrind and cover up my mistakes.
  (void) memset(dem, 0, sizeof *dem);

  dem->opts = *opts;

  dem->time.i = 0;
  dem->time.istab = 1000;
  dem->time.t = 0.0;

  dem->part.n = 0;
  dem->part.lnew = 0;

  dem->link.n = 0;

  dem->script.i = 0;
  dem->script.tprev = (double) NAN;
  dem->script.tnext = (double) NAN;

  dem->comm.i = 0;
  dem->comm.tprev = (double) NAN;
  dem->comm.tnext = (double) NAN;

  dem->cache.i = 0;
  dem->cache.tpart = (double) NAN;
  dem->cache.tprev = (double) NAN;
  dem->cache.tnext = (double) NAN;

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
      return sizeof dem->istep;
    case BMM_MSG_NUM_EKINE:
      return sizeof dem->istep + sizeof dem->est;
    case BMM_MSG_NUM_NEIGH:
      return sizeof dem->buf.neigh + sizeof dem->buf.links;
    case BMM_MSG_NUM_PARTS:
      return sizeof dem->part.n + sizeof dem->buf.parts + sizeof dem->buf.partcs;
  }

  dynamic_assert(false, "Unsupported message number");
}

bool bmm_dem_puts_stuff(struct bmm_dem const* const dem,
    enum bmm_msg_num const num) {
  switch (num) {
    case BMM_MSG_NUM_ISTEP:
      return msg_write(&dem->istep, sizeof dem->istep, NULL);
    case BMM_MSG_NUM_EKINE:
      return msg_write(&dem->istep, sizeof dem->istep, NULL) &&
        msg_write(&dem->est, sizeof dem->est, NULL);
    case BMM_MSG_NUM_NEIGH:
      return msg_write(&dem->buf.neigh, sizeof dem->buf.neigh, NULL) &&
        msg_write(&dem->buf.links, sizeof dem->buf.links, NULL);
    case BMM_MSG_NUM_PARTS:
      return msg_write(&dem->part.n, sizeof dem->part.n, NULL) &&
        msg_write(&dem->buf.parts, sizeof dem->buf.parts, NULL) &&
        msg_write(&dem->buf.partcs, sizeof dem->buf.partcs, NULL);
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

bool bmm_dem_step(struct bmm_dem* const dem) {
  switch (dem->mode) {
    case BMM_DEM_MODE_BEGIN:
      // bmm_bottom(dem);

      dem->mode = BMM_DEM_MODE_SEDIMENT;

    // TODO No! Bad touch!
    case BMM_DEM_MODE_SEDIMENT:
      if (fmod(dem->istep * dem->opts.tstep, 50.0) < dem->opts.tstep / 2.0)
        if (!bmm_disperse(dem))
          dem->mode = BMM_DEM_MODE_LINK;

      break;
    case BMM_DEM_MODE_LINK:
      if (fmod(dem->istep * dem->opts.tstep, 30.0) < dem->opts.tstep / 2.0)
        if (bmm_dem_link(dem))
          dem->mode = BMM_DEM_MODE_SMASH;

      break;
    case BMM_DEM_MODE_SMASH:
      if (fmod(dem->istep * dem->opts.tstep, 70.0) < dem->opts.tstep / 2.0)
        if (bmm_dem_break(dem)) {
          (void) bmm_dem_fix(dem);
          (void) bmm_dem_fault(dem);

          dem->mode = BMM_DEM_MODE_ACCEL;
        }

      break;
    case BMM_DEM_MODE_ACCEL:
      (void) bmm_dem_relink(dem);

      break;
  }

  if (dem->istep * dem->opts.tstep >= dem->buf.neigh.tnext) {
    bmm_dem_recont(dem);
    bmm_dem_reneigh(dem);
    // TODO Take `tstep` into account
    // since the update only happens after the drift time has passed.
    dem->buf.neigh.tnext += bmm_dem_drift(dem);
  }

  bmm_dem_predict(dem);

  bmm_dem_force(dem);

  bmm_dem_correct(dem);

  if (dem->time.i % dem->time.istab == 0)
    bmm_dem_stab(dem);

  ++dem->istep;

  return true;
}

// TODO This ought to be const-qualified.
bool bmm_dem_comm(struct bmm_dem* const dem) {
  double const tnow = dem->istep * dem->opts.tstep;

  if (tnow < dem->opts.tcomm + dem->opts.tstepcomm)
    return true;

  dem->opts.tcomm = tnow;

  // TODO This timing is bogus.
  if (dem->istep % 20 == 0)
    bmm_dem_puts(dem, BMM_MSG_NUM_NEIGH);

  // TODO Make a mechanism to automate retransmission of differences only.

  bmm_dem_puts(dem, BMM_MSG_NUM_ISTEP);
  bmm_dem_puts(dem, BMM_MSG_NUM_PARTS);

  dem->est.ekinetic = bmm_dem_ekinetic(dem);
  dem->est.pvector = bmm_dem_pvector(dem);
  dem->est.pscalar = bmm_dem_pscalar(dem);

  bmm_dem_puts(dem, BMM_MSG_NUM_EKINE);

  return true;
}

// TODO Would it be beneficial to be higher-order?
bool bmm_dem_run(struct bmm_dem* const dem) {
  int const sigs[] = {SIGINT, SIGQUIT, SIGTERM, SIGPIPE};
  if (bmm_sig_register(sigs, nmembof(sigs)) != SIZE_MAX) {
    BMM_TLE_STDS();

    return false;
  }

  bmm_dem_puts(dem, BMM_MSG_NUM_ISTEP);

#ifdef _GNU_SOURCE
#ifdef DEBUG
  int const excepts = feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  if (excepts == -1)
    BMM_TLE_STDS();
#endif
#endif

  while (dem->istep < dem->opts.nstep) {
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

    if (!bmm_dem_step(dem))
      return false;

    if (!bmm_dem_comm(dem))
      return false;
  }

#ifdef _GNU_SOURCE
#ifdef DEBUG
  if (feenableexcept(excepts) == -1)
    BMM_TLE_STDS();
#endif
#endif

  return true;
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
