#include "clopts.h"
#include "dem.h"
#include "exts.h"
#include "strs.h"
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

__attribute__ ((__nonnull__ (1, 2)))
static bool f(char const* const key, char const* const value, void* const p) {
  struct bmm_opts* const opts = p;

  if (strcmp(key, "ndim") == 0) {
    if (!bmm_strtoz(value, &opts->ndim))
      return false;
  } else if (strcmp(key, "nbin") == 0) {
    if (!bmm_strtoz(value, &opts->nbin))
      return false;
  } else if (strcmp(key, "npart") == 0) {
    if (!bmm_strtoz(value, &opts->npart))
      return false;
  } else if (strcmp(key, "nstep") == 0) {
    if (!bmm_strtoz(value, &opts->nstep))
      return false;
  }

  return true;
}

__attribute__ ((__nonnull__))
int main(int const argc, char** const argv) {
  struct bmm_opts opts;
  bmm_defopts(&opts);
  bmm_clopts((char const* const*) &argv[1], (size_t) (argc - 1), f, &opts);

  return bmm_rundem(&opts) ? EXIT_SUCCESS : EXIT_FAILURE;
}
