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
#include <math.h>
#include <signal.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
// TODO Remove this test header.
#include <unistd.h>

#ifdef _GNU_SOURCE
#ifdef DEBUG
#include <fenv.h>
#endif
#endif

extern inline void bmm_dem_clear(struct bmm_dem_list* const list);

extern inline bool bmm_dem_push(struct bmm_dem_list* const list, size_t const x);

extern inline size_t bmm_dem_size(struct bmm_dem_list const* const list);

extern inline size_t bmm_dem_get(struct bmm_dem_list const* const list, size_t const i);

extern inline struct bmm_dem_buf* bmm_dem_getbuf(struct bmm_dem*);

extern inline struct bmm_dem_buf const* bmm_dem_getrbuf(struct bmm_dem const*);

extern inline struct bmm_dem_buf* bmm_dem_getwbuf(struct bmm_dem*);

extern inline void bmm_dem_swapbuf(struct bmm_dem*);

// TODO Avoid explicitly calculating angles.

// Fake repulsive force.
void bmm_dem_fakef(struct bmm_dem* const dem) {
  struct bmm_dem_buf* const buf = bmm_dem_getbuf(dem);

  for (size_t ipart = 0; ipart < dem->opts.npart; ++ipart)
    for (size_t idim = 0; idim < 2; ++idim)
      buf->parts[ipart].lin.f[idim] = 0.0;

  /*
  for (size_t ipart = 0; ipart < dem->opts.npart; ++ipart)
    for (size_t jpart = ipart + 1; jpart < dem->opts.npart; ++jpart) {
    */
  for (size_t ipart = 0; ipart < dem->opts.npart; ++ipart)
    for (size_t ineigh = 0; ineigh < bmm_dem_size(&buf->neigh.neighs[ipart]); ++ineigh) {
      size_t const jpart = bmm_dem_get(&buf->neigh.neighs[ipart], ineigh);

      double const d2 = bmm_geom2d_pdist2(
          buf->parts[ipart].lin.r,
          buf->parts[jpart].lin.r,
          dem->rext);

      double const r2 = bmm_fp_sq(
          buf->partcs[ipart].rrad + buf->partcs[jpart].rrad);

      if (d2 < r2) {
        double const a = bmm_geom2d_pangle(
            buf->parts[ipart].lin.r,
            buf->parts[jpart].lin.r,
            dem->rext);

        double const c = 0.03;

        buf->parts[ipart].lin.f[0] -= c * cos(a);
        buf->parts[ipart].lin.f[1] -= c * sin(a);
        buf->parts[jpart].lin.f[0] += c * cos(a);
        buf->parts[jpart].lin.f[1] += c * sin(a);
      }
    }
}

void bmm_dem_euler(struct bmm_dem* const dem) {
  struct bmm_dem_buf* const buf = bmm_dem_getbuf(dem);

  // TODO Dissipate energy elsewhere.
  double const damp = 0.986;

  double const dt = dem->opts.tstep;

  for (size_t ipart = 0; ipart < dem->opts.npart; ++ipart) {
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
  struct bmm_dem_buf const* const rbuf = bmm_dem_getrbuf(dem);

  for (size_t idim = 0; idim < 2; ++idim)
    icell[idim] = bmm_size_uclamp(
        (size_t) bmm_fp_lerp(rbuf->parts[ipart].lin.r[idim],
          0.0, dem->rext[idim],
          0.0, (double) dem->opts.ncell[idim]), dem->opts.ncell[idim]);
}

void bmm_dem_recont(struct bmm_dem* const dem) {
  struct bmm_dem_buf* const buf = bmm_dem_getbuf(dem);

  for (size_t ilist = 0; ilist < bmm_size_prod(dem->opts.ncell, 2); ++ilist)
    bmm_dem_clear(&dem->pool.conts[ilist]);

  for (size_t ipart = 0; ipart < dem->opts.npart; ++ipart) {
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

  for (size_t ipart = 0; ipart < dem->opts.npart; ++ipart) {
    bmm_dem_clear(&buf->neigh.neighs[ipart]);

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
            bmm_dem_push(&buf->neigh.neighs[ipart], jpart);
        }
      }
  }
}

// Horse says neigh.
void bmm_dem_horse(struct bmm_dem* const dem) {
  struct bmm_dem_buf const* const rbuf = bmm_dem_getrbuf(dem);
  struct bmm_dem_buf* const wbuf = bmm_dem_getwbuf(dem);

  // TODO Something.
  (void) memmove(&wbuf->parts, &rbuf->parts, sizeof wbuf->parts);
  (void) memmove(&wbuf->neigh, &rbuf->neigh, sizeof wbuf->neigh);
  (void) memmove(&wbuf->partcs, &rbuf->partcs, sizeof wbuf->partcs);
}

// TODO These are dubious for empty sets.

// Maximum velocity estimator.
void bmm_dem_maxvel(double* const v, struct bmm_dem const* const dem) {
  struct bmm_dem_buf const* const rbuf = bmm_dem_getrbuf(dem);

  for (size_t idim = 0; idim < 2; ++idim) {
    v[idim] = 0.0;

    for (size_t ipart = 0; ipart < dem->opts.npart; ++ipart)
      v[idim] = fmax(v[idim], rbuf->parts[ipart].lin.v[idim]);
  }
}

// Maximum radius estimator.
double bmm_dem_maxrad(struct bmm_dem const* const dem) {
  double r = 0.0;

  struct bmm_dem_buf const* const rbuf = bmm_dem_getrbuf(dem);

  for (size_t ipart = 0; ipart < dem->opts.npart; ++ipart)
    r = fmax(r, rbuf->partcs[ipart].rrad);

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

  struct bmm_dem_buf const* const rbuf = bmm_dem_getrbuf(dem);

  for (size_t ipart = 0; ipart < dem->opts.npart; ++ipart)
    for (size_t idim = 0; idim < 2; ++idim)
      e += rbuf->partcs[ipart].mass * bmm_fp_sq(rbuf->parts[ipart].lin.v[idim]);
  // TODO No!

  return e * 0.5;
}

// Total momentum estimator.
double bmm_dem_pvector(struct bmm_dem const* const dem) {
  double p[2];
  for (size_t idim = 0; idim < 2; ++idim)
    p[idim] = 0.0;

  struct bmm_dem_buf const* const rbuf = bmm_dem_getrbuf(dem);

  for (size_t ipart = 0; ipart < dem->opts.npart; ++ipart)
    for (size_t idim = 0; idim < 2; ++idim)
      p[idim] += rbuf->partcs[ipart].mass * rbuf->parts[ipart].lin.v[idim];
  // TODO No!

  return bmm_geom2d_norm(p);
}

// Individual momentum estimator.
double bmm_dem_pscalar(struct bmm_dem const* const dem) {
  double p = 0.0;

  struct bmm_dem_buf const* const rbuf = bmm_dem_getrbuf(dem);

  for (size_t ipart = 0; ipart < dem->opts.npart; ++ipart)
    p += rbuf->partcs[ipart].mass * bmm_geom2d_norm(rbuf->parts[ipart].lin.v);
  // TODO No!

  return p;
}

void bmm_dem_defopts(struct bmm_dem_opts* const opts) {
  opts->ncell[0] = 6;
  opts->ncell[1] = 6;
  opts->nbin = 1;
  // opts->npart = 0;
  opts->npart = 256;
  opts->nstep = 2000;
  opts->rmax = 0.2;
  opts->tstep = 0.1;
  opts->vleeway = 0.01;
}

static void bmm_pretend(struct bmm_dem* const dem) {
  struct bmm_dem_buf* const buf = bmm_dem_getbuf(dem);

  for (size_t ipart = 0; ipart < dem->opts.npart; ++ipart) {
    for (size_t idim = 0; idim < 2; ++idim)
      buf->parts[ipart].lin.r[idim] = bmm_fp_uwrap(
          (double) (rand() % 256 - 128) / 128.0,
          dem->rext[idim]);

    buf->parts[ipart].ang.alpha = M_2PI * bmm_fp_uwrap(
          (double) (rand() % 256 - 128) / 128.0, 1.0);
    buf->parts[ipart].ang.omega = (double) (rand() % 256 - 128) / 128.0;
  }
}

void bmm_dem_defpartc(struct bmm_dem_partc* const partc) {
  partc->rrad = 0.0125;
  partc->mass = 1.0; // TODO No!
  partc->moi = 0.5 * partc->mass * bmm_fp_sq(partc->rrad); // TODO No!
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

void bmm_dem_def(struct bmm_dem* const dem,
    struct bmm_dem_opts const* const opts) {
  // This is here just to help Valgrind.
  memset(dem, 0, sizeof *dem);

  dem->opts = *opts;
  dem->istep = 0;

  for (size_t idim = 0; idim < 2; ++idim)
    dem->rext[idim] = 1.0;

  dem->forcesch = bmm_dem_fakef;
  dem->intsch = bmm_dem_euler;
  dem->dblbuf = false;

  if (dem->dblbuf) {
    dem->data.bufs.active = &dem->data.bufs.bufs[0];
    dem->data.bufs.passive = &dem->data.bufs.bufs[1];
  }

  struct bmm_dem_buf* const buf = bmm_dem_getbuf(dem);

  for (size_t ipart = 0; ipart < dem->opts.npart; ++ipart) {
    bmm_dem_defpartc(&buf->partcs[ipart]);
    bmm_dem_defpart(&buf->parts[ipart]);
  }

  bmm_pretend(dem); // TODO Remove later!

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
    buf->neigh.tnext += bmm_dem_drift(dem);
  }

  dem->forcesch(dem);
  dem->intsch(dem);

  if (dem->dblbuf)
    bmm_dem_horse(dem);

  bmm_dem_swapbuf(dem);

  return true;
}

static bool bmm_dem_comm(struct bmm_dem* const dem) {
  // TODO This timing is bogus.
  if (dem->istep % 20 == 0)
    bmm_putneighs(dem);

  bmm_putnop(dem);
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
  bool const result = bmm_dem_run_for_real(dem);

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
