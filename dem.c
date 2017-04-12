#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <signal.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "bit.h"
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

extern inline void bmm_dem_clearl(struct bmm_dem_listl* const list);

extern inline bool bmm_dem_pushl(struct bmm_dem_listl* const list, size_t const x);

extern inline size_t bmm_dem_sizel(struct bmm_dem_listl const* const list);

extern inline size_t bmm_dem_getl(struct bmm_dem_listl const* const list, size_t const i);

// TODO Wow, disgusting.

extern inline void bmm_dem_cleary(struct bmm_dem_listy* const list);

extern inline bool bmm_dem_pushy(struct bmm_dem_listy* const list, size_t const x);

extern inline size_t bmm_dem_sizey(struct bmm_dem_listy const* const list);

extern inline size_t bmm_dem_gety(struct bmm_dem_listy const* const list, size_t const i);

// TODO Wow, even more disgusting.

extern inline void bmm_dem_clear(struct bmm_dem_list* const list);

extern inline bool bmm_dem_push(struct bmm_dem_list* const list, size_t const x);

extern inline size_t bmm_dem_size(struct bmm_dem_list const* const list);

extern inline size_t bmm_dem_get(struct bmm_dem_list const* const list, size_t const i);

// TODO Wow, what is wrong with you?

extern inline struct bmm_dem_buf* bmm_dem_getbuf(struct bmm_dem*);

extern inline struct bmm_dem_buf const* bmm_dem_getrbuf(struct bmm_dem const*);

extern inline struct bmm_dem_buf* bmm_dem_getwbuf(struct bmm_dem*);

extern inline void bmm_dem_swapbuf(struct bmm_dem*);

// TODO Make periodicity mutable by axis.
void bmm_dem_forces(struct bmm_dem* const dem) {
  struct bmm_dem_buf* const buf = bmm_dem_getbuf(dem);

  for (size_t ipart = 0; ipart < buf->npart; ++ipart)
    for (size_t idim = 0; idim < 2; ++idim)
      buf->parts[ipart].lin.f[idim] = 0.0;

  // Cohesive (fictitious) forces.
  if (dem->mode == BMM_DEM_SEDIMENT || dem->mode == BMM_DEM_LINK)
    for (size_t ipart = 0; ipart < buf->npart; ++ipart)
      buf->parts[ipart].lin.f[1] +=
        dem->opts.fcohes * (dem->rext[1] / 2.0 - buf->parts[ipart].lin.r[1]);

  // Accelerating (currently fictitious) forces.
  if (dem->mode == BMM_DEM_ACCEL)
    for (size_t ipart = 0; ipart < buf->npart; ++ipart) {
      if (buf->parts[ipart].lin.r[1] > 0.8 * dem->rext[1])
        buf->parts[ipart].lin.f[0] += dem->opts.faccel;
      else if (buf->parts[ipart].lin.r[1] < 0.2 * dem->rext[1])
        buf->parts[ipart].lin.f[0] -= dem->opts.faccel;
    }

  // Gravitational forces.
#ifdef GRAVY
  if (dem->mode == BMM_DEM_SEDIMENT)
    for (size_t ipart = 0; ipart < buf->npart; ++ipart)
      for (size_t idim = 0; idim < 2; ++idim)
        buf->parts[ipart].lin.f[idim] +=
          dem->opts.gravy[idim] * buf->partcs[ipart].mass;
#endif

  // Pair forces.
  /*
  for (size_t ipart = 0; ipart < buf->npart; ++ipart)
    for (size_t jpart = ipart + 1; jpart < buf->npart; ++jpart) {
    */
  for (size_t ipart = 0; ipart < buf->npart; ++ipart) {
    struct bmm_dem_listy* const listy = &buf->neigh.neighs[ipart];

    for (size_t ineigh = 0; ineigh < bmm_dem_sizey(listy); ++ineigh) {
      size_t const jpart = bmm_dem_gety(listy, ineigh);

      // TODO Make a predicate function.
      double dr[2];
      bmm_geom2d_pdiff(dr,
          buf->parts[ipart].lin.r,
          buf->parts[jpart].lin.r,
          dem->rext);

      double const d2 = bmm_geom2d_norm2(dr);

      double const r = buf->partcs[ipart].rrad + buf->partcs[jpart].rrad;

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
        bmm_geom2d_rperpr(t, n);

        double df[2];
        bmm_geom2d_scale(df, n, -f);

        bmm_geom2d_add(buf->parts[ipart].lin.f, buf->parts[ipart].lin.f, df);
      }
    }
  }

  // Link forces.
  for (size_t ipart = 0; ipart < buf->npart; ++ipart) {
    struct bmm_dem_listl* const list = &buf->links[ipart];

    for (size_t ilink = 0; ilink < bmm_dem_sizel(list); ++ilink) {
      size_t const jpart = bmm_dem_getl(list, ilink);

      // TODO Copy-pasted from above...

      double dr[2];
      bmm_geom2d_pdiff(dr,
          buf->parts[ipart].lin.r,
          buf->parts[jpart].lin.r,
          dem->rext);

      double const d = bmm_geom2d_norm(dr);

      // TODO Use getters.
      double const x0 = list->linkl[ilink].x0;
      double const f = dem->opts.klink * (x0 - d);

      double n[2];
      bmm_geom2d_normal(n, dr);

      double t[2];
      bmm_geom2d_rperpr(t, n);

      double df[2];
      bmm_geom2d_scale(df, n, -f);

      bmm_geom2d_add(buf->parts[ipart].lin.f, buf->parts[ipart].lin.f, df);
    }
  }
}

void bmm_dem_euler(struct bmm_dem* const dem) {
  struct bmm_dem_buf* const buf = bmm_dem_getbuf(dem);

  double const dt = dem->opts.tstep;

  for (size_t ipart = 0; ipart < buf->npart; ++ipart)
    if (buf->partcs[ipart].free) {
      for (size_t idim = 0; idim < 2; ++idim) {
        double const a = buf->parts[ipart].lin.f[idim] / buf->partcs[ipart].mass;

        buf->parts[ipart].lin.r[idim] = bmm_fp_uwrap(
            buf->parts[ipart].lin.r[idim] +
            buf->parts[ipart].lin.v[idim] * dt, dem->rext[idim]);

        buf->parts[ipart].lin.v[idim] =
          dem->opts.damp * (buf->parts[ipart].lin.v[idim] + a * dt);
      }

      double const a = buf->parts[ipart].ang.tau / buf->partcs[ipart].moi;

      buf->parts[ipart].ang.alpha = bmm_fp_uwrap(
          buf->parts[ipart].ang.alpha +
          buf->parts[ipart].ang.omega * dt, M_2PI);

      buf->parts[ipart].ang.omega =
        dem->opts.damp * (buf->parts[ipart].ang.omega + a * dt);
    }
}

void bmm_dem_cell(size_t* const icell,
    struct bmm_dem const* const dem, size_t const ipart) {
  struct bmm_dem_buf const* const buf = bmm_dem_getrbuf(dem);

  for (size_t idim = 0; idim < 2; ++idim)
    icell[idim] = bmm_size_uclamp(
        (size_t) bmm_fp_lerp(buf->parts[ipart].lin.r[idim],
          0.0, dem->rext[idim],
          0.0, (double) dem->opts.ncell[idim]), dem->opts.ncell[idim]);
}

void bmm_dem_recont(struct bmm_dem* const dem) {
  struct bmm_dem_buf* const buf = bmm_dem_getbuf(dem);

  for (size_t ilist = 0; ilist < bmm_size_prod(dem->opts.ncell, 2); ++ilist)
    bmm_dem_clear(&dem->pool.conts[ilist]);

  for (size_t ipart = 0; ipart < buf->npart; ++ipart) {
    size_t icell[2];
    bmm_dem_cell(icell, dem, ipart);

    size_t const ilist = icell[0] + icell[1] * dem->opts.ncell[0];

    (void) bmm_dem_push(&dem->pool.conts[ilist], ipart);
  }
}

// TODO Must this always be called immediately after `bmm_dem_recont`?
// Not an important question, but interesting nevertheless.
void bmm_dem_reneigh(struct bmm_dem* const dem) {
  struct bmm_dem_buf* const buf = bmm_dem_getbuf(dem);

  for (size_t ipart = 0; ipart < buf->npart; ++ipart) {
    bmm_dem_cleary(&buf->neigh.neighs[ipart]);

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
                buf->parts[ipart].lin.r,
                buf->parts[jpart].lin.r,
                dem->rext) < bmm_fp_sq(dem->opts.rmax))
            (void) bmm_dem_pushy(&buf->neigh.neighs[ipart], jpart);
        }
      }
  }
}

// TODO Real physics.
bool bmm_dem_relink(struct bmm_dem* const dem) {
  struct bmm_dem_buf* const buf = bmm_dem_getbuf(dem);

  for (size_t ipart = 0; ipart < buf->npart; ++ipart) {
    struct bmm_dem_listl* const list = &buf->links[ipart];

    for (size_t ilink = 0; ilink < bmm_dem_sizel(list); ++ilink) {
      size_t const jpart = bmm_dem_getl(list, ilink);

      // TODO Copy-pasted from elsewhere...

      double dr[2];
      bmm_geom2d_pdiff(dr,
          buf->parts[ipart].lin.r,
          buf->parts[jpart].lin.r,
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
  struct bmm_dem_buf* const buf = bmm_dem_getbuf(dem);

  for (size_t ipart = 0; ipart < buf->npart; ++ipart) {
    struct bmm_dem_listy* const listy = &buf->neigh.neighs[ipart];

    struct bmm_dem_listl* const list = &buf->links[ipart];

    bmm_dem_clearl(list);

    for (size_t ineigh = 0; ineigh < bmm_dem_sizey(listy); ++ineigh) {
      size_t const jpart = bmm_dem_gety(listy, ineigh);

      double const d2 = bmm_geom2d_pdist2(
          buf->parts[ipart].lin.r,
          buf->parts[jpart].lin.r,
          dem->rext);

      if (d2 < bmm_fp_sq(dem->opts.linkslurp *
            (buf->partcs[ipart].rrad + buf->partcs[jpart].rrad)))
        if (bmm_dem_pushl(list, jpart))
          list->linkl[list->n - 1].x0 = sqrt(d2);
    }
  }

  return true;
}

// TODO All of it.
bool bmm_dem_break(struct bmm_dem* const dem) {
  struct bmm_dem_buf* const buf = bmm_dem_getbuf(dem);

  for (size_t ipart = 0; ipart < buf->npart; ++ipart) {
    struct bmm_dem_listl* const list = &buf->links[ipart];

    for (size_t ilink = 0; ilink < bmm_dem_sizel(list); ++ilink) {
      size_t const jpart = bmm_dem_getl(list, ilink);

      if (fabs(buf->parts[ipart].lin.r[1] - dem->rext[1] * 0.5) <
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
    if (fabs(buf->parts[ipart].lin.r[1] - dem->rext[1] * 0.5) <
        0.1 * dem->rext[1]) {
      double p = gsl_rng_uniform(dem->rng);

      if (p > 1.0) {
        // TODO Copy-pasted from up there...
        --buf->npart;
        buf->parts[ipart] = buf->parts[buf->npart];
        buf->partcs[ipart] = buf->partcs[buf->npart];
        buf->neigh.neighs[ipart] = buf->neigh.neighs[buf->npart];
        buf->links[ipart] = buf->links[buf->npart];
        --ipart;
      }
    }
  }

  return true;
}

// TODO These are dubious for empty sets.

// Maximum velocity estimator.
void bmm_dem_maxvel(double* const v, struct bmm_dem const* const dem) {
  struct bmm_dem_buf const* const buf = bmm_dem_getrbuf(dem);

  for (size_t idim = 0; idim < 2; ++idim) {
    v[idim] = 0.0;

    for (size_t ipart = 0; ipart < buf->npart; ++ipart)
      v[idim] = fmax(v[idim], buf->parts[ipart].lin.v[idim]);
  }
}

// Maximum radius estimator.
double bmm_dem_maxrad(struct bmm_dem const* const dem) {
  double r = 0.0;

  struct bmm_dem_buf const* const buf = bmm_dem_getrbuf(dem);

  for (size_t ipart = 0; ipart < buf->npart; ++ipart)
    r = fmax(r, buf->partcs[ipart].rrad);

  return r;
}

// Drift time estimator.
double bmm_dem_drift(struct bmm_dem const* const dem) {
  double t = INFINITY;

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

  struct bmm_dem_buf const* const buf = bmm_dem_getrbuf(dem);

  for (size_t ipart = 0; ipart < buf->npart; ++ipart)
    for (size_t idim = 0; idim < 2; ++idim)
      e += buf->partcs[ipart].mass * bmm_fp_sq(buf->parts[ipart].lin.v[idim]);

  return e * 0.5;
}

// Total momentum estimator.
double bmm_dem_pvector(struct bmm_dem const* const dem) {
  double p[2];
  for (size_t idim = 0; idim < 2; ++idim)
    p[idim] = 0.0;

  struct bmm_dem_buf const* const buf = bmm_dem_getrbuf(dem);

  for (size_t ipart = 0; ipart < buf->npart; ++ipart)
    for (size_t idim = 0; idim < 2; ++idim)
      p[idim] += buf->partcs[ipart].mass * buf->parts[ipart].lin.v[idim];

  return bmm_geom2d_norm(p);
}

// Individual momentum estimator.
double bmm_dem_pscalar(struct bmm_dem const* const dem) {
  double p = 0.0;

  struct bmm_dem_buf const* const buf = bmm_dem_getrbuf(dem);

  for (size_t ipart = 0; ipart < buf->npart; ++ipart)
    p += buf->partcs[ipart].mass * bmm_geom2d_norm(buf->parts[ipart].lin.v);

  return p;
}

// Mean coefficient of restitution
// (just linear dashpot for now, also a bit wrong).
double bmm_dem_cor(struct bmm_dem const* const dem) {
  double e = 0.0;

  struct bmm_dem_buf const* const buf = bmm_dem_getrbuf(dem);

  for (size_t ipart = 0; ipart < buf->npart; ++ipart) {
    double const meff =
      buf->partcs[ipart].mass * buf->partcs[ipart].mass /
      (buf->partcs[ipart].mass + buf->partcs[ipart].mass);
    e += exp(-M_PI * dem->opts.yelast / (2.0 * meff) /
        sqrt(dem->opts.ymodul / meff -
          bmm_fp_sq(dem->opts.yelast / (2.0 * meff))));
  }

  return e / (double) buf->npart;
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
  opts->klink = 120.0;
  opts->fcohes = 0.02;
  opts->faccel = 0.01;
  opts->gravy[0] = 0.0;
  opts->gravy[1] = -0.005;
  opts->ymodul = 2.0e+4;
  opts->yelast = 1.4e+2;
  // TODO Dissipate energy elsewhere.
  opts->damp = 0.999;
  opts->rmean = 0.04;
  opts->lucky = 13;

#ifdef NDEBUG
  for (size_t idim = 0; idim < 2; ++idim)
    opts->ncell[idim] *= 2;
  opts->tstep /= 2.0;
  opts->klink *= 4.0;
  opts->fcohes /= 2.0;
  opts->rmean /= 2.0;
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
  struct bmm_dem_buf* const buf = bmm_dem_getbuf(dem);

  double x = 0.0;

  for ever {
    double r = dem->opts.rmean + gsl_ran_gaussian(dem->rng, dem->opts.rsd);

    x += r;

    struct bmm_dem_partc* const partc = &buf->partcs[buf->npart];
    struct bmm_dem_part* const part = &buf->parts[buf->npart];

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

      ++buf->npart;

      break;
    }

    x += r;

    ++buf->npart;
  }
}

static bool bmm_disperse(struct bmm_dem* const dem) {
  size_t success = 0;

  struct bmm_dem_buf* const buf = bmm_dem_getbuf(dem);

  size_t maxfail = 100;
  size_t fail = 0;

  while (fail < maxfail && buf->npart < BMM_PART_MAX) {
    double const r = dem->opts.rmean +
      gsl_ran_gaussian(dem->rng, dem->opts.rsd);

    double x[2];
    x[0] = gsl_rng_uniform(dem->rng) * dem->rext[0];
    x[1] = gsl_rng_uniform(dem->rng) * dem->rext[1];

    for (size_t ipart = 0; ipart < buf->npart; ++ipart) {
      double dr[2];
      bmm_geom2d_pdiff(dr,
          buf->parts[ipart].lin.r,
          x,
          dem->rext);

      double const d2 = bmm_geom2d_norm2(dr);

      double const rs = buf->partcs[ipart].rrad + r;

      double const r2 = bmm_fp_sq(rs);

      if (d2 < r2) {
        ++fail;
        goto ret;
      }
    }

    fail = 0;

    struct bmm_dem_partc* const partc = &buf->partcs[buf->npart];
    struct bmm_dem_part* const part = &buf->parts[buf->npart];

    bmm_dem_partc_def(partc);
    bmm_dem_part_def(part);

    partc->rrad = r;
    partc->free = true;
    part->lin.r[0] = x[0];
    part->lin.r[1] = x[1];

    ++buf->npart;

    ++success;
ret: ;
  }

  // TODO The lucky number.
  return success > dem->opts.lucky;
}

void bmm_dem_def(struct bmm_dem* const dem,
    struct bmm_dem_opts const* const opts) {
  // This is here just to help Valgrind and cover up my mistakes.
  memset(dem, 0, sizeof *dem);

  dem->opts = *opts;
  dem->mode = BMM_DEM_BEGIN;
  dem->istep = 0;

  for (size_t idim = 0; idim < 2; ++idim)
    dem->rext[idim] = 1.0;

  dem->forcesch = bmm_dem_forces;
  dem->intsch = bmm_dem_euler;
  dem->dblbuf = false;

  if (dem->dblbuf) {
    dem->data.bufs.active = &dem->data.bufs.bufs[0];
    dem->data.bufs.passive = &dem->data.bufs.bufs[1];
  }

  struct bmm_dem_buf* const buf = bmm_dem_getbuf(dem);

  buf->npart = 0;

  buf->neigh.tnext = 0.0;
}

// TODO Is this error handling bad?
static bool msg_read(uint8_t const* buf, size_t const n,
    __attribute__ ((__unused__)) void* const ptr) {
  switch (bmm_io_readin(buf, n)) {
    case BMM_IO_READ_ERROR:
    case BMM_IO_READ_EOF:
      return false;
    case BMM_IO_READ_SUCCESS:
      return true;
  }
}

static bool msg_write(uint8_t const* buf, size_t const n,
    __attribute__ ((__unused__)) void* const ptr) {
  return bmm_io_writeout(buf, n);
}

static size_t bmm_dem_sniff_size(struct bmm_dem const* const dem,
    enum bmm_msg_type const type) {
  struct bmm_dem_buf const* const buf = bmm_dem_getrbuf(dem);

  switch (type) {
    case BMM_MSG_NPART:
      return sizeof buf->npart;
    case BMM_MSG_EKINE:
      return sizeof dem->istep + sizeof dem->est;
  }

  dynamic_assert(false, "Unsupported message type");
}

static bool bmm_dem_gets_stuff(struct bmm_dem* const dem,
    enum bmm_msg_type const type, size_t const size) {
  struct bmm_dem_buf* const buf = bmm_dem_getbuf(dem);

  switch (type) {
    case BMM_MSG_NPART:
      return bmm_msg_data_read(&buf->npart, msg_read, size);
    case BMM_MSG_EKINE:
      return bmm_msg_data_read(&dem->istep, msg_read, size) &&
        bmm_msg_data_write(&dem->est, msg_write, size);
  }

  dynamic_assert(false, "Unsupported message type");
}

static bool bmm_dem_gets(struct bmm_dem* const dem,
    enum bmm_msg_type* const type) {
  struct bmm_msg_spec spec;
  if (!bmm_msg_spec_read(&spec, msg_read, NULL))
    return false;

  if (spec.endian != BMM_MSG_ENDIAN_LITTLE) {
    BMM_TLE_EXTS(BMM_TLE_UNIMPL, "Unsupported endianness");

    return false;
  }

  if (spec.tag != BMM_MSG_TAG_SP) {
    BMM_TLE_EXTS(BMM_TLE_UNIMPL, "Unsupported tag");

    return false;
  }

  size_t const asize = spec.msg.size - BMM_MSG_TYPESIZE;

  if (!bmm_msg_type_read(type, msg_read, NULL))
    return false;

  size_t const esize = bmm_dem_sniff_size(dem, *type);

  if (esize != asize)
    BMM_TLE_EXTS(BMM_TLE_UNKNOWN, "Size mismatch");
  else if (esize < asize) {
    BMM_TLE_EXTS(BMM_TLE_UNKNOWN, "Buffer would overflow");

    return false;
  }

  return bmm_dem_gets_stuff(dem, *type, asize);
}

static bool bmm_dem_puts_stuff(struct bmm_dem const* const dem,
    enum bmm_msg_type const type, size_t const size) {
  struct bmm_dem_buf const* const buf = bmm_dem_getrbuf(dem);

  switch (type) {
    case BMM_MSG_NPART:
      return bmm_msg_data_write(&buf->npart, msg_write, size);
    case BMM_MSG_EKINE:
      return bmm_msg_data_write(&dem->istep, msg_write, size) &&
        bmm_msg_data_write(&dem->est, msg_write, size);
  }

  dynamic_assert(false, "Unsupported message type");
}

static bool bmm_dem_puts(struct bmm_dem const* const dem,
    enum bmm_msg_type const type) {
  size_t const size = bmm_dem_sniff_size(dem, type);

  struct bmm_msg_spec spec;
  bmm_msg_spec_def(&spec);
  spec.endian = BMM_MSG_ENDIAN_LITTLE;
  spec.msg.size = size + BMM_MSG_TYPESIZE;

  return bmm_msg_spec_write(&spec, msg_write, NULL) &&
    bmm_msg_type_write(&type, msg_write, NULL) &&
    bmm_dem_puts_stuff(dem, type, size);
}

static void bmm_dem_put(struct bmm_dem const* const dem,
    enum bmm_msg_type const msg) {
  struct bmm_msg_head head;
  bmm_head_def(&head);
  bmm_bit_pset(&head.flags, BMM_FBIT_INTLE);
  bmm_bit_pset(&head.flags, BMM_FBIT_FPLE);
  bmm_bit_pset(&head.flags, BMM_FBIT_FLUSH);

  head.type = msg;

  bmm_msg_put(&head, dem);
}

bool bmm_dem_step(struct bmm_dem* const dem) {
  struct bmm_dem_buf* const buf = bmm_dem_getbuf(dem);

  switch (dem->mode) {
    case BMM_DEM_BEGIN:
      // bmm_bottom(dem);

      dem->mode = BMM_DEM_SEDIMENT;

    // TODO No! Bad touch!
    case BMM_DEM_SEDIMENT:
      if (fmod(dem->istep * dem->opts.tstep, 50.0) < dem->opts.tstep / 2.0)
        if (!bmm_disperse(dem))
          dem->mode = BMM_DEM_LINK;

      break;
    case BMM_DEM_LINK:
      if (fmod(dem->istep * dem->opts.tstep, 30.0) < dem->opts.tstep / 2.0)
        if (bmm_dem_link(dem))
          dem->mode = BMM_DEM_BREAK;

      break;
    case BMM_DEM_BREAK:
      if (fmod(dem->istep * dem->opts.tstep, 70.0) < dem->opts.tstep / 2.0)
        if (bmm_dem_break(dem))
          dem->mode = BMM_DEM_ACCEL;

      break;
    case BMM_DEM_ACCEL:
      (void) bmm_dem_relink(dem);

      break;
  }

  if (dem->istep * dem->opts.tstep >= buf->neigh.tnext) {
    bmm_dem_recont(dem);
    bmm_dem_reneigh(dem);
    buf->neigh.tnext += bmm_dem_drift(dem);
  }

  dem->forcesch(dem);
  dem->intsch(dem);

  bmm_dem_swapbuf(dem);

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
    bmm_dem_put(dem, BMM_MSG_NEIGH);

  // TODO Make a mechanism to automate retransmission of differences only.

  bmm_dem_put(dem, BMM_MSG_NPART);
  bmm_dem_put(dem, BMM_MSG_PARTS);

  dem->est.ekinetic = bmm_dem_ekinetic(dem);
  dem->est.pvector = bmm_dem_pvector(dem);
  dem->est.pscalar = bmm_dem_pscalar(dem);

  bmm_dem_put(dem, BMM_MSG_EKINE);

  return true;
}

// TODO Would it be beneficial to be higher-order?
bool bmm_dem_run(struct bmm_dem* const dem) {
  int const sigs[] = {SIGINT, SIGQUIT, SIGTERM, SIGPIPE};
  if (bmm_sig_register(sigs, sizeof sigs / sizeof *sigs) != SIZE_MAX) {
    BMM_TLE_STDS();

    return false;
  }

  bmm_dem_put(dem, BMM_MSG_NPART);

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
          BMM_TLE_EXTS(BMM_TLE_ASYNC, "Simulation interrupted");

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

static bool bmm_dem_run_rng(struct bmm_dem* const dem) {
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

  bool const result = bmm_dem_run_rng(dem);

  free(dem);

  return result;
}
