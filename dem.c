#include "bit.h"
#include "conf.h"
#include "dem.h"
#include "err.h"
#include "fp.h"
#include "geom.h"
#include "geom2d.h"
#include "msg.h"
#include "sig.h"
#include "size.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <signal.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifdef _GNU_SOURCE
#ifdef DEBUG
#include <fenv.h>
#endif
#endif

extern inline void bmm_dem_cleary(struct bmm_dem_listy* const list);

extern inline bool bmm_dem_pushy(struct bmm_dem_listy* const list, size_t const x);

extern inline size_t bmm_dem_sizey(struct bmm_dem_listy const* const list);

extern inline size_t bmm_dem_gety(struct bmm_dem_listy const* const list, size_t const i);

// TODO Wow, disgusting.

extern inline void bmm_dem_clear(struct bmm_dem_list* const list);

extern inline bool bmm_dem_push(struct bmm_dem_list* const list, size_t const x);

extern inline size_t bmm_dem_size(struct bmm_dem_list const* const list);

extern inline size_t bmm_dem_get(struct bmm_dem_list const* const list, size_t const i);

// TODO Wow, even more disgusting.

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

  if (dem->mode == BMM_DEM_SEDIMENT)
    for (size_t ipart = 0; ipart < buf->npart; ++ipart)
      for (size_t idim = 0; idim < 2; ++idim)
        buf->parts[ipart].lin.f[idim] +=
          dem->opts.gravy[idim] * buf->partcs[ipart].mass;

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
        double const dx = listy->thingy[ineigh].x[1] - listy->thingy[ineigh].x[0];
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
}

void bmm_dem_euler(struct bmm_dem* const dem) {
  struct bmm_dem_buf* const buf = bmm_dem_getbuf(dem);

  // TODO Dissipate energy elsewhere.
  double const damp = 0.986;

  double const dt = dem->opts.tstep;

  for (size_t ipart = 0; ipart < buf->npart; ++ipart)
    if (buf->partcs[ipart].free) {
      for (size_t idim = 0; idim < 2; ++idim) {
        double const a = buf->parts[ipart].lin.f[idim] / buf->partcs[ipart].mass;

        buf->parts[ipart].lin.r[idim] = bmm_fp_uwrap(
            buf->parts[ipart].lin.r[idim] +
            buf->parts[ipart].lin.v[idim] * dt, dem->rext[idim]);

        buf->parts[ipart].lin.v[idim] =
          damp * (buf->parts[ipart].lin.v[idim] + a * dt);
      }

      double const a = buf->parts[ipart].ang.tau / buf->partcs[ipart].moi;

      buf->parts[ipart].ang.alpha = bmm_fp_uwrap(
          buf->parts[ipart].ang.alpha +
          buf->parts[ipart].ang.omega * dt, M_2PI);

      buf->parts[ipart].ang.omega =
        damp * (buf->parts[ipart].ang.omega + a * dt);
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

void bmm_dem_relink(struct bmm_dem* const dem) {
  struct bmm_dem_buf* const buf = bmm_dem_getbuf(dem);

  for (size_t ipart = 0; ipart < buf->npart; ++ipart) {
    struct bmm_dem_listy* const listy = &buf->neigh.neighs[ipart];

    struct bmm_dem_list* const list = &buf->links[ipart];

    bmm_dem_clear(list);

    for (size_t ineigh = 0; ineigh < bmm_dem_sizey(listy); ++ineigh) {
      size_t const jpart = bmm_dem_gety(listy, ineigh);

      // TODO This is bogus; use triangulation instead.
      if (bmm_geom2d_pdist2(
            buf->parts[ipart].lin.r,
            buf->parts[jpart].lin.r,
            dem->rext) < bmm_fp_sq(2.0 *
              (dem->opts.rmean + dem->opts.rsd)))
        (void) bmm_dem_push(list, jpart);
    }
  }
}

// Horse says neigh.
void bmm_dem_horse(struct bmm_dem* const dem) {
  struct bmm_dem_buf const* const buf = bmm_dem_getrbuf(dem);
  struct bmm_dem_buf* const wbuf = bmm_dem_getwbuf(dem);

  // TODO Something.
  (void) memmove(&wbuf->parts, &buf->parts, sizeof wbuf->parts);
  (void) memmove(&wbuf->neigh, &buf->neigh, sizeof wbuf->neigh);
  (void) memmove(&wbuf->partcs, &buf->partcs, sizeof wbuf->partcs);
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

void bmm_dem_defopts(struct bmm_dem_opts* const opts) {
  opts->ncell[0] = 6;
  opts->ncell[1] = 6;
  opts->nbin = 1;
  opts->nstep = 100000;
  opts->rmax = 0.2;
  // opts->tend = ...;
  // opts->tadv = ...;
  opts->tstep = 0.01;
  opts->tstepcomm = 1.0;
  opts->tcomm = 0.0;
  opts->vleeway = 0.01;
  opts->gravy[0] = 0.0;
  opts->gravy[1] = -0.005;
  opts->ymodul = 2.0e+4;
  opts->yelast = 1.0e+2;
  opts->rmean = 0.03;
  opts->rsd = opts->rmean * 0.2;
}

void bmm_dem_defpartc(struct bmm_dem_partc* const partc) {
  partc->rrad = 0.0125;
  partc->mass = 1.0;
  partc->moi = bmm_geom_ballmoi(partc->rrad, 3) * partc->mass;
  partc->free = true;
}

void bmm_dem_defpart(struct bmm_dem_part* const part) {
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

    bmm_dem_defpartc(partc);
    bmm_dem_defpart(part);

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

static void bmm_disperse(struct bmm_dem* const dem) {
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
        goto end;
      }
    }

    fail = 0;

    struct bmm_dem_partc* const partc = &buf->partcs[buf->npart];
    struct bmm_dem_part* const part = &buf->parts[buf->npart];

    bmm_dem_defpartc(partc);
    bmm_dem_defpart(part);

    partc->rrad = r;
    partc->free = true;
    part->lin.r[0] = x[0];
    part->lin.r[1] = x[1];

    ++buf->npart;
end: ;
  }
}

void bmm_dem_def(struct bmm_dem* const dem,
    struct bmm_dem_opts const* const opts) {
  // This is here just to help Valgrind and cover up my mistakes.
  memset(dem, 0, sizeof *dem);

  dem->opts = *opts;
  dem->mode = BMM_DEM_SEDIMENT;
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

// TODO Unify these three.

static void bmm_putnop(struct bmm_dem const* const dem) {
  struct bmm_msg_head head;
  bmm_defhead(&head);
  bmm_bit_pset(&head.flags, BMM_FBIT_INTLE);
  bmm_bit_pset(&head.flags, BMM_FBIT_FPLE);

  head.type = BMM_MSG_NOP;

  bmm_msg_put(&head, dem);
}

static void bmm_putopts(struct bmm_dem const* const dem) {
  struct bmm_msg_head head;
  bmm_defhead(&head);
  bmm_bit_pset(&head.flags, BMM_FBIT_INTLE);
  bmm_bit_pset(&head.flags, BMM_FBIT_FPLE);

  head.type = BMM_MSG_NPART;

  bmm_msg_put(&head, dem);
}

static void bmm_putparts(struct bmm_dem const* const dem) {
  struct bmm_msg_head head;
  bmm_defhead(&head);
  bmm_bit_pset(&head.flags, BMM_FBIT_INTLE);
  bmm_bit_pset(&head.flags, BMM_FBIT_FPLE);
  bmm_bit_pset(&head.flags, BMM_FBIT_FLUSH);

  head.type = BMM_MSG_PARTS;

  bmm_msg_put(&head, dem);
}

static void bmm_putneighs(struct bmm_dem const* const dem) {
  struct bmm_msg_head head;
  bmm_defhead(&head);
  bmm_bit_pset(&head.flags, BMM_FBIT_INTLE);
  bmm_bit_pset(&head.flags, BMM_FBIT_FPLE);
  bmm_bit_pset(&head.flags, BMM_FBIT_FLUSH);

  head.type = BMM_MSG_NEIGH;

  bmm_msg_put(&head, dem);
}

static bool bmm_dem_step(struct bmm_dem* const dem) {
  struct bmm_dem_buf* const buf = bmm_dem_getbuf(dem);

  if (dem->istep * dem->opts.tstep >= buf->neigh.tnext) {
    bmm_dem_recont(dem);
    bmm_dem_reneigh(dem);
    bmm_dem_relink(dem);
    buf->neigh.tnext += bmm_dem_drift(dem);
  }

  // TODO Make a mechanism to automate retransmission of differences.

  dem->forcesch(dem);
  dem->intsch(dem);

  if (dem->dblbuf)
    bmm_dem_horse(dem);

  bmm_dem_swapbuf(dem);

  return true;
}

static bool bmm_dem_comm(struct bmm_dem* const dem) {
  double const tnow = dem->istep * dem->opts.tstep;

  if (tnow < dem->opts.tcomm + dem->opts.tstepcomm)
    return true;

  dem->opts.tcomm = tnow;

  // TODO This timing is bogus.
  if (dem->istep % 20 == 0)
    bmm_putneighs(dem);

  bmm_putnop(dem);
  bmm_putopts(dem);
  bmm_putparts(dem);

  // TODO These should go via messages.
  dem->est.ekinetic = bmm_dem_ekinetic(dem);
  dem->est.pvector = bmm_dem_pvector(dem);
  dem->est.pscalar = bmm_dem_pscalar(dem);
  // fprintf(stderr, "%f\n", bmm_dem_ekinetic(dem));
  // fprintf(stderr, "%f\n", bmm_dem_pvector(dem));

  return true;
}

static bool bmm_dem_run_for_real(struct bmm_dem* const dem) {
  bmm_putopts(dem);

  // TODO Remove these test messages.
  while (dem->istep < dem->opts.nstep) {
    int signum;
    if (bmm_sig_use(&signum))
      switch (signum) {
        case SIGINT:
        case SIGQUIT:
        case SIGTERM:
        case SIGPIPE:
          BMM_ERR_FWARN(NULL, "Simulation interrupted");

          return false;
      }

    if (!bmm_dem_step(dem))
      return false;

    if (!bmm_dem_comm(dem))
      return false;

    ++dem->istep;
  }

  return true;
}

static bool bmm_dem_run_now(struct bmm_dem_opts const* const opts) {
  struct bmm_dem* const dem = malloc(sizeof *dem);
  if (dem == NULL) {
    BMM_ERR_WARN(malloc);

    return false;
  }

  bmm_dem_def(dem, opts);

  // TODO Stack frame.
  gsl_rng_type const* const t = gsl_rng_env_setup();
  dem->rng = gsl_rng_alloc(t);
  if (dem->rng == NULL) {
    BMM_ERR_WARN(gsl_rng_alloc);

    free(dem);
    return false;
  }

  bmm_bottom(dem);

  bmm_disperse(dem);

  bool const result = bmm_dem_run_for_real(dem);

  gsl_rng_free(dem->rng);

  free(dem);

  return result;
}

bool bmm_dem_run(struct bmm_dem_opts const* const opts) {
  int const sigs[] = {SIGINT, SIGQUIT, SIGTERM, SIGPIPE};
  if (bmm_sig_register(sigs, sizeof sigs / sizeof *sigs) != SIZE_MAX) {
    BMM_ERR_WARN(bmm_sig_register);

    return false;
  }

#ifdef _GNU_SOURCE
#ifdef DEBUG
  int const excepts = feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  if (excepts == -1)
    BMM_ERR_WARN(feenableexcept);
#endif
#endif

  bool const result = bmm_dem_run_now(opts);

#ifdef _GNU_SOURCE
#ifdef DEBUG
  if (feenableexcept(excepts) == -1)
    BMM_ERR_WARN(feenableexcept);
#endif
#endif

  return result;
}
