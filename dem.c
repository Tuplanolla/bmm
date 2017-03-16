#include "bit.h"
#include "conf.h"
#include "dem.h"
#include "err.h"
#include "msg.h"
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
// TODO Remove this test header.
#include <unistd.h>

#ifdef _GNU_SOURCE
#ifdef DEBUG
#include <fenv.h>
#endif
#endif

// TODO Unify these two.

void bmm_putopts(struct bmm_dem const* const dem) {
  struct bmm_head head;
  bmm_defhead(&head);
  bmm_bit_pset(&head.flags, BMM_FBIT_INTLE);
  bmm_bit_pset(&head.flags, BMM_FBIT_FPLE);

  head.type = BMM_MSG_NDIM;

  bmm_msg_put(&head, dem);
}

void bmm_putparts(struct bmm_dem const* const dem) {
  struct bmm_head head;
  bmm_defhead(&head);
  bmm_bit_pset(&head.flags, BMM_FBIT_INTLE);
  bmm_bit_pset(&head.flags, BMM_FBIT_FPLE);
  bmm_bit_pset(&head.flags, BMM_FBIT_FLUSH);

  head.type = BMM_MSG_PARTS;

  bmm_msg_put(&head, dem);
}

void bmm_pretend(struct bmm_dem* const dem) {
  for (size_t ipart = 0; ipart < BMM_PART_MAX; ++ipart)
    for (size_t idim = 0; idim < BMM_DIM_MAX; ++idim)
      dem->parts[ipart].rpos[idim] += (double) (rand() % 256 - 128) * 100e-6;
}

void bmm_dem_defopts(struct bmm_dem_opts* const opts) {
  opts->ndim = 0;
  opts->nbin = 0;
  opts->npart = 0;
  opts->nstep = 0;
}

void bmm_dem_defpart(struct bmm_dem_part* const part) {
  part->rrad = 0.0;
  part->arot = 0.0;

  for (size_t idim = 0; idim < BMM_DIM_MAX; ++idim)
    part->rpos[idim] = 0.0;
}

void bmm_dem_def(struct bmm_dem* const dem,
    struct bmm_dem_opts const* const opts) {
  dem->opts = *opts;
  dem->istep = 0;

  for (size_t idim = 0; idim < BMM_DIM_MAX; ++idim)
    dem->rexts[idim] = 1.0;

  for (size_t ipart = 0; ipart < BMM_PART_MAX; ++ipart)
    bmm_dem_defpart(&dem->parts[ipart]);
}

bool bmm_dem_run(struct bmm_dem_opts const* const opts) {
  struct bmm_dem dem;
  bmm_dem_def(&dem, opts);
  dem.opts = *opts;

#ifdef _GNU_SOURCE
#ifdef DEBUG
  int const excepts = feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  if (excepts == -1)
    BMM_ERR_WARN(feenableexcept);
#endif
#endif

  bmm_putopts(&dem);

  // TODO Remove these test messages.
  for (size_t istep = 0; istep < dem.opts.nstep; ++istep) {
    bmm_pretend(&dem);
    bmm_putparts(&dem);
    // usleep(100000);
  }

#ifdef _GNU_SOURCE
#ifdef DEBUG
  if (feenableexcept(excepts) == -1)
    BMM_ERR_WARN(feenableexcept);
#endif
#endif

  return true;
}