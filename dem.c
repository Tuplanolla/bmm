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

bool bmm_rundem(struct bmm_opts const* const opts) {
  struct bmm_state state;
  state.opts = *opts;

#ifdef _GNU_SOURCE
#ifdef DEBUG
  int const excepts = feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  if (excepts == -1)
    bmm_error("Failed to enable floating-point exceptions.");
#endif
#endif

  // TODO Remove this test message.
  struct bmm_head head;
  head.type = BMM_MSG_NDIM;
  bmm_putmsg(stdout, &head, &state);

#ifdef _GNU_SOURCE
#ifdef DEBUG
  if (feenableexcept(excepts) == -1)
    bmm_error("Failed to restore floating-point exceptions.");
#endif
#endif

  return true;
}
