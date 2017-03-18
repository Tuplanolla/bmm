#include "bit.h"
#include "conf.h"
#include "dem.h"
#include "err.h"
#include "msg.h"
#include "sig.h"
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

extern inline struct bmm_dem_buf const* bmm_dem_getrbuf(struct bmm_dem const*);

extern inline struct bmm_dem_buf* bmm_dem_getwbuf(struct bmm_dem*);

extern inline void bmm_dem_swapbuf(struct bmm_dem*);

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

static void bmm_pretend(struct bmm_dem* const dem) {
  struct bmm_dem_buf* const buf = bmm_dem_getwbuf(dem);

  for (size_t ipart = 0; ipart < BMM_PART_MAX; ++ipart)
    for (size_t idim = 0; idim < 2; ++idim)
      buf->parts[ipart].rpos[idim] += (double) (rand() % 256 - 128) * 100e-6;

  bmm_dem_swapbuf(dem);
}

void bmm_dem_defopts(struct bmm_dem_opts* const opts) {
  opts->ncell[0] = 1;
  opts->ncell[1] = 1;
  opts->nbin = 0;
  opts->npart = 0;
  opts->nstep = 0;
}

void bmm_dem_defpart(struct bmm_dem_part* const part) {
  part->rrad = 0.0;
  part->arot = 0.0;

  for (size_t idim = 0; idim < 2; ++idim)
    part->rpos[idim] = 0.0;
}

void bmm_dem_def(struct bmm_dem* const dem,
    struct bmm_dem_opts const* const opts) {
  dem->opts = *opts;
  dem->istep = 0;

  for (size_t idim = 0; idim < 2; ++idim)
    dem->rexts[idim] = 1.0;

  dem->dblbuf = false;

  struct bmm_dem_buf* const buf = bmm_dem_getwbuf(dem);

  for (size_t ipart = 0; ipart < BMM_PART_MAX; ++ipart)
    bmm_dem_defpart(&buf->parts[ipart]);
}

static bool bmm_dem_run_for_real(struct bmm_dem* const dem) {
  bmm_putopts(dem);

  // TODO Remove these test messages.
  for (size_t istep = 0; istep < dem->opts.nstep; ++istep) {
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

    bmm_putnop(dem);
    bmm_pretend(dem);
    bmm_putparts(dem);
    // usleep(100000);
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
