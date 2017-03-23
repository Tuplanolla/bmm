#include "bit.h"
#include "conf.h"
#include "dem.h"
#include "err.h"
#include "fp.h"
#include "msg.h"
#include "sig.h"
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

  for (size_t ipart = 0; ipart < dem->opts.npart; ++ipart)
    for (size_t jpart = ipart + 1; jpart < dem->opts.npart; ++jpart) {
      double const d2 = bmm_fp_pdist2(
          buf->parts[ipart].lin.r,
          buf->parts[jpart].lin.r,
          dem->rext);

      double const r2 = bmm_fp_sq(
          buf->parts[ipart].rrad + buf->parts[jpart].rrad);

      if (d2 < r2) {
        double const a = bmm_fp_pangle(
            buf->parts[ipart].lin.r,
            buf->parts[jpart].lin.r,
            dem->rext);

        double const c = 0.04;

        buf->parts[ipart].lin.f[0] -= c * cos(a);
        buf->parts[ipart].lin.f[1] -= c * sin(a);
        buf->parts[jpart].lin.f[0] += c * cos(a);
        buf->parts[jpart].lin.f[1] += c * sin(a);
      }
    }
}

void bmm_dem_euler(struct bmm_dem* const dem) {
  struct bmm_dem_buf const* const rbuf = bmm_dem_getrbuf(dem);
  struct bmm_dem_buf* const wbuf = bmm_dem_getwbuf(dem);

  double const dt = dem->tstep;

  for (size_t ipart = 0; ipart < dem->opts.npart; ++ipart) {
    // TODO No!
    wbuf->parts[ipart].rrad = rbuf->parts[ipart].rrad;
    wbuf->parts[ipart].mass = rbuf->parts[ipart].mass;

    for (size_t idim = 0; idim < 2; ++idim) {
      double const a = rbuf->parts[ipart].lin.f[idim] / rbuf->parts[ipart].mass;

      wbuf->parts[ipart].lin.v[idim] =
        rbuf->parts[ipart].lin.v[idim] + a * dt;

      wbuf->parts[ipart].lin.r[idim] = bmm_fp_uwrap(
          rbuf->parts[ipart].lin.r[idim] +
          rbuf->parts[ipart].lin.v[idim] * dt, dem->rext[idim]);
    }

    // TODO This is bogus (needs to be the moment of inertia).
    double const a = rbuf->parts[ipart].ang.tau / rbuf->parts[ipart].mass;

    wbuf->parts[ipart].ang.omega =
      rbuf->parts[ipart].ang.omega + a * dt;

    wbuf->parts[ipart].ang.alpha = bmm_fp_uwrap(
        rbuf->parts[ipart].ang.alpha +
        rbuf->parts[ipart].ang.omega * dt, M_2PI);
  }
}

// Horse says neigh.
void bmm_dem_horse(struct bmm_dem* const dem) {
  struct bmm_dem_buf const* const rbuf = bmm_dem_getrbuf(dem);
  struct bmm_dem_buf* const wbuf = bmm_dem_getwbuf(dem);

  // TODO Something.
  (void) memmove(&wbuf->neigh, &rbuf->neigh, sizeof wbuf->neigh);
  (void) memmove(&wbuf->partcs, &rbuf->partcs, sizeof wbuf->partcs);
}

// Total kinetic energy estimator.
double bmm_dem_ekinetic(struct bmm_dem const* const dem) {
  double e = 0.0;

  struct bmm_dem_buf const* const rbuf = bmm_dem_getrbuf(dem);

  for (size_t ipart = 0; ipart < dem->opts.npart; ++ipart)
    for (size_t idim = 0; idim < 2; ++idim)
      e += rbuf->parts[ipart].mass * bmm_fp_sq(rbuf->parts[ipart].lin.v[idim]);
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
      p[idim] += rbuf->parts[ipart].mass * rbuf->parts[ipart].lin.v[idim];
  // TODO No!

  return bmm_fp_norm(p);
}

// Individual momentum estimator.
double bmm_dem_pscalar(struct bmm_dem const* const dem) {
  double p = 0.0;

  struct bmm_dem_buf const* const rbuf = bmm_dem_getrbuf(dem);

  for (size_t ipart = 0; ipart < dem->opts.npart; ++ipart)
    p += rbuf->parts[ipart].mass * bmm_fp_norm(rbuf->parts[ipart].lin.v);
  // TODO No!

  return p;
}

void bmm_dem_defopts(struct bmm_dem_opts* const opts) {
  opts->ncell[0] = 4;
  opts->ncell[1] = 4;
  opts->nbin = 1;
  // opts->npart = 0;
  opts->npart = 256;
  opts->nstep = 2000;
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

void bmm_dem_defpart(struct bmm_dem_part* const part) {
  part->rrad = 0.0125;
  part->mass = 1.0; // TODO No!

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
  dem->tstep = 0.1;

  for (size_t idim = 0; idim < 2; ++idim)
    dem->rext[idim] = 1.0;

  dem->forcesch = bmm_dem_fakef;
  dem->intsch = bmm_dem_euler;
  dem->dblbuf = true;

  if (dem->dblbuf) {
    dem->data.bufs.active = &dem->data.bufs.bufs[0];
    dem->data.bufs.passive = &dem->data.bufs.bufs[1];
  }

  struct bmm_dem_buf* const buf = bmm_dem_getbuf(dem);

  for (size_t ipart = 0; ipart < dem->opts.npart; ++ipart)
    bmm_dem_defpart(&buf->parts[ipart]);

  bmm_pretend(dem); // TODO Remove later!
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

static bool bmm_dem_step(struct bmm_dem* const dem) {
  dem->forcesch(dem);
  dem->intsch(dem);
  bmm_dem_horse(dem);

  bmm_dem_swapbuf(dem);

  return true;
}

static bool bmm_dem_comm(struct bmm_dem* const dem) {
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
