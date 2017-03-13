#include "bits.h"
#include "dem.h"
#include "errors.h"
#include "msgs.h"
#include <stdbool.h>
#include <stdlib.h>

#ifdef _GNU_SOURCE
#ifdef DEBUG
#include <fenv.h>
#endif
#endif

// TODO Unify these two.

void bmm_putopts(struct bmm_state const* const state) {
  struct bmm_head head;
  bmm_defhead(&head);
  bmm_setbp(&head.flags, BMM_FBIT_INTLE);
  bmm_setbp(&head.flags, BMM_FBIT_FPLE);

  head.type = BMM_MSG_NDIM;

  bmm_putmsg(stdout, &head, state);
}

void bmm_putparts(struct bmm_state const* const state) {
  struct bmm_head head;
  bmm_defhead(&head);
  bmm_setbp(&head.flags, BMM_FBIT_INTLE);
  bmm_setbp(&head.flags, BMM_FBIT_FPLE);

  head.type = BMM_MSG_PARTS;

  bmm_putmsg(stdout, &head, state);
}

void bmm_defopts(struct bmm_opts* const opts) {
  opts->ndim = 0;
  opts->nbin = 0;
  opts->npart = 0;
  opts->nstep = 0;
}

void bmm_defpart(struct bmm_part* const part) {
  part->rrad = 0.0;

  for (size_t idim = 0; idim < DIM_MAX; ++idim)
    part->rpos[idim] = 0.0;
}

void bmm_defstate(struct bmm_state* const state) {
  bmm_defopts(&state->opts);
  state->istep = 0;

  for (size_t idim = 0; idim < DIM_MAX; ++idim)
    state->rexts[idim] = 0.0;

  for (size_t ipart = 0; ipart < PART_MAX; ++ipart)
    bmm_defpart(&state->parts[ipart]);
}

bool bmm_rundem(struct bmm_opts const* const opts) {
  struct bmm_state state;
  bmm_defstate(&state);
  state.opts = *opts;

#ifdef _GNU_SOURCE
#ifdef DEBUG
  int const excepts = feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  if (excepts == -1)
    bmm_error("Failed to enable floating-point exceptions.");
#endif
#endif

  // TODO Remove these test messages.
  bmm_putopts(&state);
  bmm_putparts(&state);

#ifdef _GNU_SOURCE
#ifdef DEBUG
  if (feenableexcept(excepts) == -1)
    bmm_error("Failed to restore floating-point exceptions.");
#endif
#endif

  return true;
}
