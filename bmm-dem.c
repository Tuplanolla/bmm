#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "dem.h"
#include "ext.h"
#include "opt.h"
#include "str.h"
#include "tle.h"

__attribute__ ((__nonnull__ (1, 2)))
static bool f(char const* const key, char const* const value,
    void* const ptr) {
  struct bmm_dem_opts* const opts = ptr;

  if (strcmp(key, "ncellx") == 0) {
    if (!bmm_str_strtoz(&opts->ncell[0], value))
      return false;
  } else if (strcmp(key, "ncelly") == 0) {
    if (!bmm_str_strtoz(&opts->ncell[1], value))
      return false;
  } else if (strcmp(key, "nbin") == 0) {
    if (!bmm_str_strtoz(&opts->nbin, value))
      return false;
  } else if (strcmp(key, "nstep") == 0) {
    if (!bmm_str_strtoz(&opts->nstep, value))
      return false;
  } else if (strcmp(key, "rmax") == 0) {
    if (!bmm_str_strtod(&opts->rmax, value))
      return false;
  } else
    return false;

  return true;
}

__attribute__ ((__nonnull__))
int main(int const argc, char** const argv) {
  bmm_tle_reset(argv[0]);

  struct bmm_dem_opts opts;
  bmm_dem_opts_def(&opts);

  if (!bmm_opt_parse((char const* const*) &argv[1], (size_t) (argc - 1),
        f, &opts)) {
    bmm_tle_put();

    return EXIT_FAILURE;
  }

  // TODO Get rid of these after refactoring `dem.c`.
  opts->cache.ncell[0] = 6;
  opts->cache.ncell[1] = 6;
  opts->cache.rcutoff = 0.2;

  if (!bmm_dem_run_with(&opts)) {
    bmm_tle_put();

    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
