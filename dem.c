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
#include "msg.h"
#include "sig.h"
#include "size.h"
#include "tle.h"

#ifdef _GNU_SOURCE
#ifdef DEBUG
#include <fenv.h>
#endif
#endif

// The partition looks like this:
//
//       ^
//       | - + +
//     y | - + +
//       | - - +
//       +------->
//           x
//
// The following monad comprehension generator was used:
//
//     [(i, j) | i <- [-1 .. 1], j <- [-1 .. 1]]
int const bmm_moore[][2] = {
  {-1, -1}, {-1, 0}, {-1, 1},
  {0, -1}, {0, 0}, {0, 1},
  {1, -1}, {1, 0}, {1, 1}
};

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

    dynamic_assert(n >= 3, "Too few neighbor cells");

    size_t const j = n - 1;
    double const r = (double) j - 1.0;
    double const x = bmm_fp_lerp(dem->part.r[ipart][idim],
        0.0, dem->opts.box.r[idim], 0.0, r);

    if (x < 0.0)
      pijcell[idim] = 0;
    else if (x >= r)
      pijcell[idim] = j;
    else {
      size_t const k = (size_t) x;

      dynamic_assert(k < j, "Incorrect rounding");

      pijcell[idim] = k;
    }
  }
}

// TODO Test iso and specialize.

void bmm_dem_i2ijcell(size_t* const pijcell,
    struct bmm_dem const* const dem, size_t icell) {
  for (size_t idim = 0; idim < BMM_NDIM; ++idim) {
    bmm_size_div_t const qr = bmm_size_div(icell,
        dem->opts.cache.ncell[BMM_NDIM - 1 - idim]);

    icell = qr.quot;
    pijcell[idim] = qr.rem;
  }
}

size_t bmm_dem_ij2icell(struct bmm_dem* const dem,
    size_t const* const ijcell) {
  size_t icell = 0;

  for (size_t idim = 0; idim < BMM_NDIM; ++idim) {
    icell *= dem->opts.cache.ncell[BMM_NDIM - 1 - idim];
    icell += ijcell[idim];
  }

  return icell;
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

  dem->part.phi[ipart][idim] = 0.0;
  dem->part.omega[ipart][idim] = 0.0;
  dem->part.alpha[ipart][idim] = 0.0;

  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    dem->part.f[ipart][idim] = 0.0;

  dem->part.tau[ipart][idim] = 0.0;

  // TODO Calculate local neighbors, because this is as fast as a bulk update.
  //
  // Usually you fill `cache.part` by going over all particles.
  // Then you go over all particles and,
  // for each particle in the same cell or one Moore half,
  // check the neighborhood condition.
  //
  // In this case you only go over this particle and,
  // for each particle in the same cell or one Moore half,
  // check the neighborhood condition.
  // However you also go through every particle in the other Moore half,
  // checking the neighborhood only for this particle.

  // "In this case..."
  for (size_t imoore = (BMM_POW(3, BMM_NDIM) - 1) / 2;
      imoore < BMM_POW(3, BMM_NDIM);
      ++imoore) {
  }

  // TODO Would like this interface.
  size_t const n = bmm_moore_nge(BMM_NDIM,
      dem->opts.box.per, dem->opts.cache.ncell);

  for (size_t imoore = 0; imoore < n; ++imoore) {
    size_t const icell = bmm_moore_icell(imoore, BMM_NDIM,
        dem->opts.box.per, dem->opts.cache.ncell);

    foreach particle in icell
      if (bmm_geom2d_cpdist2(
            dem->buf.parts[ipart].lin.r,
            dem->buf.parts[jpart].lin.r,
            dem->opts.box.per, dem->opts.box.r) < bmm_fp_sq(dem->opts.rmax))
        (void) bmm_dem_pushy(&dem->buf.neigh.neighs[ipart], jpart);
  }

  // "However..."
  for (size_t imoore = 0;
      imoore < (BMM_POW(3, BMM_NDIM) - 1) / 2;
      ++imoore) {
  }

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
    size_t const icell = bmm_dem_ij2icell(dem, dem->cache.cell[jpart]);

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

// TODO Make periodicity mutable by axis.
void bmm_dem_forces(struct bmm_dem* const dem) {
  for (size_t ipart = 0; ipart < dem->buf.npart; ++ipart)
    for (size_t idim = 0; idim < 2; ++idim)
      dem->buf.parts[ipart].lin.f[idim] = 0.0;

  // Cohesive (fictitious) forces.
  if (dem->mode == BMM_DEM_MODE_SEDIMENT || dem->mode == BMM_DEM_MODE_LINK)
    for (size_t ipart = 0; ipart < dem->buf.npart; ++ipart)
      dem->buf.parts[ipart].lin.f[1] +=
        dem->opts.fcohes * (dem->rext[1] / 2.0 - dem->buf.parts[ipart].lin.r[1]);

  // Gravitational forces.
#ifdef GRAVY
  if (dem->mode == BMM_DEM_MODE_SEDIMENT)
    for (size_t ipart = 0; ipart < dem->buf.npart; ++ipart)
      for (size_t idim = 0; idim < 2; ++idim)
        dem->buf.parts[ipart].lin.f[idim] +=
          dem->opts.gravy[idim] * dem->buf.partcs[ipart].mass;
#endif

  // Pair forces.
  /*
  for (size_t ipart = 0; ipart < dem->buf.npart; ++ipart)
    for (size_t jpart = ipart + 1; jpart < dem->buf.npart; ++jpart) {
    */
  for (size_t ipart = 0; ipart < dem->buf.npart; ++ipart) {
    struct bmm_dem_listy* const listy = &dem->buf.neigh.neighs[ipart];

    for (size_t ineigh = 0; ineigh < bmm_dem_sizey(listy); ++ineigh) {
      size_t const jpart = bmm_dem_gety(listy, ineigh);

      // TODO Make a predicate function.
      double dr[2];
      bmm_geom2d_pdiff(dr,
          dem->buf.parts[ipart].lin.r,
          dem->buf.parts[jpart].lin.r,
          dem->rext);

      double const d2 = bmm_geom2d_norm2(dr);

      double const r = dem->buf.partcs[ipart].rrad + dem->buf.partcs[jpart].rrad;

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
            dem->opts.ymodul * x + dem->opts.yelast * v);

        double n[2];
        bmm_geom2d_normal(n, dr);

        double t[2];
        bmm_geom2d_rperp(t, n);

        double df[2];
        bmm_geom2d_scale(df, n, -f);

        bmm_geom2d_add(dem->buf.parts[ipart].lin.f, dem->buf.parts[ipart].lin.f, df);
      }
    }
  }

  // Link forces.
  for (size_t ipart = 0; ipart < dem->buf.npart; ++ipart) {
    struct bmm_dem_listl* const list = &dem->buf.links[ipart];

    for (size_t ilink = 0; ilink < bmm_dem_sizel(list); ++ilink) {
      size_t const jpart = bmm_dem_getl(list, ilink);

      // TODO Copy-pasted from above...

      double dr[2];
      bmm_geom2d_pdiff(dr,
          dem->buf.parts[ipart].lin.r,
          dem->buf.parts[jpart].lin.r,
          dem->rext);

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

      bmm_geom2d_add(dem->buf.parts[ipart].lin.f, dem->buf.parts[ipart].lin.f, df);
    }
  }

  // TODO Calculate force feedback from the residuals.

  // Accelerating (currently fictitious) forces.
  if (dem->mode == BMM_DEM_MODE_ACCEL) {
    // TODO Make an estimator for this driven total velocity?

    size_t ntotal = 0;
    double vtotal[2];
    for (size_t idim = 0; idim < 2; ++idim)
      vtotal[idim] = 0.0;

    for (size_t ipart = 0; ipart < dem->buf.npart; ++ipart)
      if (!dem->buf.partcs[ipart].nondr) {
        ++ntotal;

        for (size_t idim = 0; idim < 2; ++idim)
          vtotal[idim] += dem->buf.parts[ipart].lin.v[idim];
      }

    for (size_t idim = 0; idim < 2; ++idim) {
      if (vtotal[idim] / (double) ntotal < dem->opts.vcrunch[idim])
        dem->fcrunch[idim] += dem->opts.fadjust;
      else
        dem->fcrunch[idim] -= dem->opts.fadjust;
    }

    for (size_t ipart = 0; ipart < dem->buf.npart; ++ipart) {
      if (!dem->buf.partcs[ipart].free)
        for (size_t idim = 0; idim < 2; ++idim)
          dem->buf.parts[ipart].lin.f[idim] = 0.0;
      if (!dem->buf.partcs[ipart].nondr)
        for (size_t idim = 0; idim < 2; ++idim)
          dem->buf.parts[ipart].lin.f[idim] += dem->fcrunch[idim];
    }
  }
}

void bmm_dem_euler(struct bmm_dem* const dem) {
  double const dt = dem->opts.tstep;

  for (size_t ipart = 0; ipart < dem->buf.npart; ++ipart) {
    for (size_t idim = 0; idim < 2; ++idim) {
      double const a = dem->buf.parts[ipart].lin.f[idim] / dem->buf.partcs[ipart].mass;

      dem->buf.parts[ipart].lin.r[idim] = bmm_fp_uwrap(
          dem->buf.parts[ipart].lin.r[idim] +
          dem->buf.parts[ipart].lin.v[idim] * dt, dem->rext[idim]);

      dem->buf.parts[ipart].lin.v[idim] =
        dem->opts.damp * (dem->buf.parts[ipart].lin.v[idim] + a * dt);
    }

    double const a = dem->buf.parts[ipart].ang.tau / dem->buf.partcs[ipart].moi;

    dem->buf.parts[ipart].ang.alpha = bmm_fp_uwrap(
        dem->buf.parts[ipart].ang.alpha +
        dem->buf.parts[ipart].ang.omega * dt, M_2PI);

    dem->buf.parts[ipart].ang.omega =
      dem->opts.damp * (dem->buf.parts[ipart].ang.omega + a * dt);
  }
}

void bmm_dem_cell(size_t* const icell,
    struct bmm_dem const* const dem, size_t const ipart) {
  for (size_t idim = 0; idim < 2; ++idim)
    icell[idim] = bmm_size_uclamp(
        (size_t) bmm_fp_lerp(dem->buf.parts[ipart].lin.r[idim],
          0.0, dem->rext[idim],
          0.0, (double) dem->opts.ncell[idim]), dem->opts.ncell[idim]);
}

void bmm_dem_recont(struct bmm_dem* const dem) {
  for (size_t ilist = 0; ilist < bmm_size_prod(dem->opts.ncell, 2); ++ilist)
    bmm_dem_clear(&dem->pool.conts[ilist]);

  for (size_t ipart = 0; ipart < dem->buf.npart; ++ipart) {
    size_t icell[2];
    bmm_dem_cell(icell, dem, ipart);

    size_t const ilist = icell[0] + icell[1] * dem->opts.ncell[0];

    (void) bmm_dem_push(&dem->pool.conts[ilist], ipart);
  }
}

// TODO Must this always be called immediately after `bmm_dem_recont`?
// Not an important question, but interesting nevertheless.
void bmm_dem_reneigh(struct bmm_dem* const dem) {
  for (size_t ipart = 0; ipart < dem->buf.npart; ++ipart) {
    bmm_dem_cleary(&dem->buf.neigh.neighs[ipart]);

    size_t icell[2];
    bmm_dem_cell(icell, dem, ipart);

    size_t (* const f[])(size_t, size_t) = {
      bmm_size_dec, bmm_size_constant, bmm_size_inc
    };

    // The size is $3^d$.
    for (size_t imoore = 0; imoore < 3; ++imoore)
      for (size_t jmoore = 0; jmoore < 3; ++jmoore) {
        size_t const iwhat = f[imoore](icell[0], dem->opts.ncell[0]);
        size_t const jwhat = f[jmoore](icell[1], dem->opts.ncell[1]);

        size_t const icont = iwhat + jwhat * dem->opts.ncell[0];

        for (size_t ineigh = 0; ineigh < bmm_dem_size(&dem->pool.conts[icont]); ++ineigh) {
          size_t const jpart = bmm_dem_get(&dem->pool.conts[icont], ineigh);

          if (bmm_geom2d_pdist2(
                dem->buf.parts[ipart].lin.r,
                dem->buf.parts[jpart].lin.r,
                dem->rext) < bmm_fp_sq(dem->opts.rmax))
            (void) bmm_dem_pushy(&dem->buf.neigh.neighs[ipart], jpart);
        }
      }
  }
}

// TODO Real physics.
bool bmm_dem_relink(struct bmm_dem* const dem) {
  for (size_t ipart = 0; ipart < dem->buf.npart; ++ipart) {
    struct bmm_dem_listl* const list = &dem->buf.links[ipart];

    for (size_t ilink = 0; ilink < bmm_dem_sizel(list); ++ilink) {
      size_t const jpart = bmm_dem_getl(list, ilink);

      // TODO Copy-pasted from elsewhere...

      double dr[2];
      bmm_geom2d_pdiff(dr,
          dem->buf.parts[ipart].lin.r,
          dem->buf.parts[jpart].lin.r,
          dem->rext);

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
  for (size_t ipart = 0; ipart < dem->buf.npart; ++ipart) {
    struct bmm_dem_listy* const listy = &dem->buf.neigh.neighs[ipart];

    struct bmm_dem_listl* const list = &dem->buf.links[ipart];

    bmm_dem_clearl(list);

    for (size_t ineigh = 0; ineigh < bmm_dem_sizey(listy); ++ineigh) {
      size_t const jpart = bmm_dem_gety(listy, ineigh);

      if (jpart == ipart)
        continue;

      double const d2 = bmm_geom2d_pdist2(
          dem->buf.parts[ipart].lin.r,
          dem->buf.parts[jpart].lin.r,
          dem->rext);

      if (d2 < bmm_fp_sq(dem->opts.linkslurp *
            (dem->buf.partcs[ipart].rrad + dem->buf.partcs[jpart].rrad)))
        if (bmm_dem_pushl(list, jpart))
          list->linkl[list->n - 1].x0 = dem->opts.linkoff * sqrt(d2);
    }
  }

  return true;
}

// TODO All of it.
bool bmm_dem_break(struct bmm_dem* const dem) {
  for (size_t ipart = 0; ipart < dem->buf.npart; ++ipart) {
    struct bmm_dem_listl* const list = &dem->buf.links[ipart];

    for (size_t ilink = 0; ilink < bmm_dem_sizel(list); ++ilink) {
      size_t const jpart = bmm_dem_getl(list, ilink);

      if (fabs(dem->buf.parts[ipart].lin.r[1] - dem->rext[1] * 0.5) <
          0.1 * dem->rext[1]) {
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
    if (fabs(dem->buf.parts[ipart].lin.r[1] - dem->rext[1] * 0.5) <
        0.1 * dem->rext[1]) {
      double p = gsl_rng_uniform(dem->rng);

      if (p > 1.0) {
        // TODO Copy-pasted from up there...
        --dem->buf.npart;
        dem->buf.parts[ipart] = dem->buf.parts[dem->buf.npart];
        dem->buf.partcs[ipart] = dem->buf.partcs[dem->buf.npart];
        dem->buf.neigh.neighs[ipart] = dem->buf.neigh.neighs[dem->buf.npart];
        dem->buf.links[ipart] = dem->buf.links[dem->buf.npart];
        --ipart;
      }
    }
  }

  return true;
}

bool bmm_dem_fix(struct bmm_dem* const dem) {
  for (size_t ipart = 0; ipart < dem->buf.npart; ++ipart) {
    if (dem->buf.parts[ipart].lin.r[1] > 0.8 * dem->rext[1]) {
      dem->buf.partcs[ipart].free = false;
      dem->buf.partcs[ipart].nondr = false;
    } else if (dem->buf.parts[ipart].lin.r[1] < 0.2 * dem->rext[1])
      dem->buf.partcs[ipart].free = false;
  }

  return true;
}

bool bmm_dem_fault(struct bmm_dem* const dem) {
  for (size_t ipart = 0; ipart < dem->buf.npart; ++ipart) {
    if (dem->buf.parts[ipart].lin.r[1] > 0.5 * dem->rext[1])
      dem->buf.parts[ipart].lin.r[1] += dem->opts.yoink * dem->rext[1];
    else
      dem->buf.parts[ipart].lin.r[1] -= dem->opts.yoink * dem->rext[1];
  }

  return true;
}

// TODO These are dubious for empty sets.

// Maximum velocity estimator.
void bmm_dem_maxvel(double* const v, struct bmm_dem const* const dem) {
  for (size_t idim = 0; idim < 2; ++idim) {
    v[idim] = 0.0;

    for (size_t ipart = 0; ipart < dem->buf.npart; ++ipart)
      v[idim] = fmax(v[idim], dem->buf.parts[ipart].lin.v[idim]);
  }
}

// Maximum radius estimator.
double bmm_dem_maxrad(struct bmm_dem const* const dem) {
  double r = 0.0;

  for (size_t ipart = 0; ipart < dem->buf.npart; ++ipart)
    r = fmax(r, dem->buf.partcs[ipart].rrad);

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
            0.5 * dem->rext[idim] / (double) dem->opts.ncell[idim]) - rad) /
        (v[idim] + dem->opts.vleeway));

  return t;
}

// Total kinetic energy estimator.
double bmm_dem_ekinetic(struct bmm_dem const* const dem) {
  double e = 0.0;

  for (size_t ipart = 0; ipart < dem->buf.npart; ++ipart)
    for (size_t idim = 0; idim < 2; ++idim)
      e += dem->buf.partcs[ipart].mass * bmm_fp_sq(dem->buf.parts[ipart].lin.v[idim]);

  return e * 0.5;
}

// Total momentum estimator.
double bmm_dem_pvector(struct bmm_dem const* const dem) {
  double p[2];
  for (size_t idim = 0; idim < 2; ++idim)
    p[idim] = 0.0;

  for (size_t ipart = 0; ipart < dem->buf.npart; ++ipart)
    for (size_t idim = 0; idim < 2; ++idim)
      p[idim] += dem->buf.partcs[ipart].mass * dem->buf.parts[ipart].lin.v[idim];

  return bmm_geom2d_norm(p);
}

// Individual momentum estimator.
double bmm_dem_pscalar(struct bmm_dem const* const dem) {
  double p = 0.0;

  for (size_t ipart = 0; ipart < dem->buf.npart; ++ipart)
    p += dem->buf.partcs[ipart].mass * bmm_geom2d_norm(dem->buf.parts[ipart].lin.v);

  return p;
}

// Mean coefficient of restitution
// (just linear dashpot for now, also a bit wrong).
double bmm_dem_cor(struct bmm_dem const* const dem) {
  double e = 0.0;

  for (size_t ipart = 0; ipart < dem->buf.npart; ++ipart) {
    double const meff =
      dem->buf.partcs[ipart].mass * dem->buf.partcs[ipart].mass /
      (dem->buf.partcs[ipart].mass + dem->buf.partcs[ipart].mass);
    e += exp(-M_PI * dem->opts.yelast / (2.0 * meff) /
        sqrt(dem->opts.ymodul / meff -
          bmm_fp_sq(dem->opts.yelast / (2.0 * meff))));
  }

  return e / (double) dem->buf.npart;
}

void bmm_dem_opts_def(struct bmm_dem_opts* const opts) {
  opts->ncell[0] = 6;
  opts->ncell[1] = 6;
  opts->nbin = 1;
  opts->rmax = 0.2;
  // opts->tend = ...;
  // opts->tadv = ...;
  opts->tstep = 0.004;
  opts->tstepcomm = 1.0;
  opts->tcomm = 0.0;
  opts->vleeway = 0.01;
  opts->linkslurp = 1.1;
  opts->linkoff = 0.95;
  opts->klink = 120.0;
  opts->fcohes = 0.02;
  opts->fadjust = 1e-7;
  opts->vcrunch[0] = 0.005;
  opts->vcrunch[1] = 0.0;
  opts->gravy[0] = 0.0;
  opts->gravy[1] = -0.005;
  opts->ymodul = 2.0e+4;
  opts->yelast = 1.4e+2;
  // TODO Dissipate energy elsewhere.
  opts->damp = 0.999;
  opts->rmean = 0.04;
  opts->yoink = 0.01;
  opts->lucky = 13;

#ifdef NDEBUG
  for (size_t idim = 0; idim < 2; ++idim)
    opts->ncell[idim] *= 2;
  opts->tstep /= 2.0;
  opts->klink *= 4.0;
  opts->fcohes /= 2.0;
  opts->rmean /= 2.0;
  opts->vcrunch[0] /= 2.0;
  opts->lucky *= 8;
#endif

  opts->nstep = (size_t) (500 / opts->tstep); // ??
  opts->rsd = opts->rmean * 0.2;
}

void bmm_dem_partc_def(struct bmm_dem_partc* const partc) {
  partc->rrad = 0.0125;
  partc->mass = 1.0;
  partc->moi = bmm_geom_ballmoi(partc->rrad, 3) * partc->mass;
  partc->free = true;
  partc->nondr = true;
}

void bmm_dem_part_def(struct bmm_dem_part* const part) {
  for (size_t idim = 0; idim < 2; ++idim) {
    part->lin.r[idim] = 0.5;
    part->lin.v[idim] = 0.0;
    part->lin.f[idim] = 0.0;
  }

  part->ang.alpha = 0.0;
  part->ang.omega = 0.0;
  part->ang.tau = 0.0;
}

// TODO Rewrite these.

static void bmm_bottom(struct bmm_dem* const dem) {
  double x = 0.0;

  for ever {
    double r = dem->opts.rmean + gsl_ran_gaussian(dem->rng, dem->opts.rsd);

    x += r;

    struct bmm_dem_partc* const partc = &dem->buf.partcs[dem->buf.npart];
    struct bmm_dem_part* const part = &dem->buf.parts[dem->buf.npart];

    bmm_dem_partc_def(partc);
    bmm_dem_part_def(part);

    partc->rrad = r;
    partc->free = false;
    part->lin.r[0] = x;
    part->lin.r[1] = dem->rext[1] / 32.0;

    if (x + r >= dem->rext[0]) {
      x -= r;

      double const rprime = (dem->rext[0] - x) / 2.0;

      x += rprime;

      partc->rrad = rprime;
      partc->free = false;
      part->lin.r[0] = x;
      part->lin.r[1] = dem->rext[1] / 32.0;

      ++dem->buf.npart;

      break;
    }

    x += r;

    ++dem->buf.npart;
  }
}

static bool bmm_disperse(struct bmm_dem* const dem) {
  size_t success = 0;

  size_t maxfail = 100;
  size_t fail = 0;

  while (fail < maxfail && dem->buf.npart < BMM_NPART) {
    double const r = dem->opts.rmean +
      gsl_ran_gaussian(dem->rng, dem->opts.rsd);

    double x[2];
    x[0] = gsl_rng_uniform(dem->rng) * dem->rext[0];
    x[1] = gsl_rng_uniform(dem->rng) * dem->rext[1];

    for (size_t ipart = 0; ipart < dem->buf.npart; ++ipart) {
      double dr[2];
      bmm_geom2d_pdiff(dr,
          dem->buf.parts[ipart].lin.r,
          x,
          dem->rext);

      double const d2 = bmm_geom2d_norm2(dr);

      double const rs = dem->buf.partcs[ipart].rrad + r;

      double const r2 = bmm_fp_sq(rs);

      if (d2 < r2) {
        ++fail;
        goto ret;
      }
    }

    fail = 0;

    struct bmm_dem_partc* const partc = &dem->buf.partcs[dem->buf.npart];
    struct bmm_dem_part* const part = &dem->buf.parts[dem->buf.npart];

    bmm_dem_partc_def(partc);
    bmm_dem_part_def(part);

    partc->rrad = r;
    partc->free = true;
    part->lin.r[0] = x[0];
    part->lin.r[1] = x[1];

    ++dem->buf.npart;

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
  dem->mode = BMM_DEM_MODE_BEGIN;
  dem->istep = 0;

  for (size_t idim = 0; idim < 2; ++idim)
    dem->rext[idim] = 1.0;

  dem->forcesch = bmm_dem_forces;
  dem->intsch = bmm_dem_euler;

  for (size_t idim = 0; idim < 2; ++idim)
    dem->fcrunch[idim] = 0.0;

  dem->buf.npart = 0;

  dem->buf.neigh.tnext = 0.0;
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
      return sizeof dem->buf.npart + sizeof dem->buf.parts + sizeof dem->buf.partcs;
  }

  dynamic_assert(false, "Unsupported message num");
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
      return msg_write(&dem->buf.npart, sizeof dem->buf.npart, NULL) &&
        msg_write(&dem->buf.parts, sizeof dem->buf.parts, NULL) &&
        msg_write(&dem->buf.partcs, sizeof dem->buf.partcs, NULL);
  }

  dynamic_assert(false, "Unsupported message num");
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
          dem->mode = BMM_DEM_MODE_BREAK;

      break;
    case BMM_DEM_MODE_BREAK:
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

  dem->forcesch(dem);
  dem->intsch(dem);

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
