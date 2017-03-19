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

// TODO Move over to another translation unit.

static void bmm_dem_diff(double* const rdiff,
    double const* const r0, double const* const r1,
    double const* const rexts) {
  for (size_t idim = 0; idim < 2; ++idim)
    rdiff[idim] = bmm_fp_swrap(r1[idim] - r0[idim], rexts[idim]);
}

// This maps $(x, y)$ to $(r^2, \phi)$.
// Goes fast.
inline void bmm_fp_to_polar2(double* const polar2, double const* const cart) {
  double r2 = 0.0;
  for (size_t idim = 0; idim < 2; ++idim)
    r2 += bmm_fp_sq(cart[idim]);

  polar2[0] = r2;
  polar2[1] = atan2(cart[1], cart[0]);
}

// This maps $(x, y)$ to $(r, \phi)$.
// Goes slow.
inline void bmm_fp_to_polar(double* const polar, double const* const cart) {
  bmm_fp_to_polar2(polar, cart);
  polar[0] = sqrt(polar[0]);
}

// This maps $(r, \phi)$ to $(x, y)$.
// Goes fast.
inline void bmm_fp_from_polar(double* const cart, double const* const polar) {
  cart[0] = polar[0] * cos(polar[1]);
  cart[1] = polar[0] * sin(polar[1]);
}

// This maps $(r^2, \phi)$ to $(x, y)$.
// Goes slow.
inline void bmm_fp_from_polar2(double* const cart,
    double const* const polar2) {
  double const polar[] = {sqrt(polar2[0]), polar2[1]};
  bmm_fp_from_polar(cart, polar);
}

// TODO No shit.

// Shit distance.
static double bmm_dem_dist2(double const* const r0, double const* const r1) {
  double d = 0.0;
  for (size_t idim = 0; idim < 2; ++idim)
    d += bmm_fp_sq(r1[idim] - r0[idim]);

  return d;
}

// Shit angle.
static double bmm_dem_angle(double const* const r0, double const* const r1) {
  double dr[2];
  for (size_t idim = 0; idim < 2; ++idim)
    dr[idim] = r1[idim] - r0[idim];

  return atan2(dr[1], dr[0]);
}

// Shit periodic distance (minimum image convention).
static double bmm_dem_pdist2(double const* const rexts,
    double const* const r0, double const* const r1) {
  double d = 0.0;

  for (size_t idim = 0; idim < 2; ++idim)
    d += bmm_fp_sq(bmm_fp_swrap(r1[idim] - r0[idim], rexts[idim]));

  return d;
}

void bmm_dem_fakef(struct bmm_dem* const dem) {
  struct bmm_dem_buf* const buf = bmm_dem_getbuf(dem);

  for (size_t ipart = 0; ipart < dem->opts.npart; ++ipart)
    for (size_t jpart = ipart + 1; jpart < dem->opts.npart; ++jpart) {
      double const d2 = bmm_dem_dist2(
          buf->parts[ipart].lin.r,
          buf->parts[jpart].lin.r);

      if (d2 < bmm_fp_sq(buf->parts[ipart].rrad + buf->parts[jpart].rrad)) {
        double const a = bmm_dem_angle(
            buf->parts[ipart].lin.r,
            buf->parts[jpart].lin.r);

        buf->parts[ipart].lin.f[0] -= d2 * cos(a);
        buf->parts[ipart].lin.f[1] -= d2 * sin(a);
        buf->parts[jpart].lin.f[0] += d2 * cos(a);
        buf->parts[jpart].lin.f[1] += d2 * sin(a);
      }
    }
}

// TODO Buffering woes.
void bmm_dem_euler(struct bmm_dem* const dem) {
  struct bmm_dem_buf const* const rbuf = bmm_dem_getrbuf(dem);
  struct bmm_dem_buf* const wbuf = bmm_dem_getwbuf(dem);

  double const dt = dem->tstep;

  for (size_t ipart = 0; ipart < dem->opts.npart; ++ipart)
    for (size_t idim = 0; idim < 2; ++idim) {
      double const a = rbuf->parts[ipart].lin.f[idim] / rbuf->parts[ipart].mass;

      wbuf->parts[ipart].lin.r[idim] += rbuf->parts[ipart].lin.v[idim] * dt;
      wbuf->parts[ipart].lin.v[idim] += a * dt;
    }

  bmm_dem_swapbuf(dem);
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

void bmm_dem_defopts(struct bmm_dem_opts* const opts) {
  opts->ncell[0] = 1;
  opts->ncell[1] = 1;
  opts->nbin = 1;
  // opts->npart = 0;
  opts->npart = 8;
  opts->nstep = 60;
}

static void bmm_pretend(struct bmm_dem* const dem) {
  struct bmm_dem_buf* const buf = bmm_dem_getbuf(dem);

  for (size_t ipart = 0; ipart < dem->opts.npart; ++ipart)
    for (size_t idim = 0; idim < 2; ++idim)
      buf->parts[ipart].lin.r[idim] += (double) (rand() % 256 - 128) * 1e-3;
}

void bmm_dem_defpart(struct bmm_dem_part* const part) {
  part->rrad = 0.1;
  part->mass = 1.0; // TODO No!

  part->ang.alpha = 0.0;
  part->ang.omega = 0.0;
  part->ang.tau = 0.0;

  for (size_t idim = 0; idim < 2; ++idim) {
    part->lin.r[idim] = 0.0;
    part->lin.v[idim] = 0.0;
    part->lin.f[idim] = 0.0;
  }
}

void bmm_dem_def(struct bmm_dem* const dem,
    struct bmm_dem_opts const* const opts) {
  dem->opts = *opts;
  dem->istep = 0;
  dem->tstep = 0.01;

  for (size_t idim = 0; idim < 2; ++idim)
    dem->rexts[idim] = 1.0;

  dem->forcesch = bmm_dem_fakef;
  dem->intsch = bmm_dem_euler;
  dem->dblbuf = false;

  if (dem->dblbuf) {
    dem->data.bufs.active = &dem->data.bufs.bufs[0];
    dem->data.bufs.passive = &dem->data.bufs.bufs[1];
  }

  struct bmm_dem_buf* const buf = bmm_dem_getbuf(dem);

  for (size_t ipart = 0; ipart < dem->opts.npart; ++ipart)
    bmm_dem_defpart(&buf->parts[ipart]);

  bmm_pretend(dem); // TODO Remove later!
}

static bool bmm_dem_step(struct bmm_dem* const dem) {
  dem->forcesch(dem);
  dem->intsch(dem);

  return true;
}

static bool bmm_dem_comm(struct bmm_dem* const dem) {
  bmm_putnop(dem);
  bmm_putparts(dem);

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
