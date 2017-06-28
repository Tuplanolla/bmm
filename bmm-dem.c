#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#ifdef _GNU_SOURCE
#ifdef DEBUG
#include <fenv.h>
#endif
#endif

#include "dem.h"
#include "ext.h"
#include "opt.h"
#include "str.h"
#include "tle.h"

static double const dtstuff = 1.0e-4;

__attribute__ ((__deprecated__))
static bool bmm_dem_script_pushidle(struct bmm_dem_opts* const opts,
    double const tspan) {
  if (opts->script.n >= nmembof(opts->script.mode))
    return false;

  opts->script.tspan[opts->script.n] = tspan;
  opts->script.dt[opts->script.n] = dtstuff;
  opts->script.mode[opts->script.n] = BMM_DEM_MODE_IDLE;

  ++opts->script.n;

  return true;
}

__attribute__ ((__deprecated__))
static bool bmm_dem_script_pushcreate(struct bmm_dem_opts* const opts,
    double const eta) {
  if (opts->script.n >= nmembof(opts->script.mode))
    return false;

  opts->script.tspan[opts->script.n] = 0.0;
  opts->script.dt[opts->script.n] = 1.0e-9;
  opts->script.mode[opts->script.n] = BMM_DEM_MODE_CREATE;
  opts->script.params[opts->script.n].create.eta = eta;

  ++opts->script.n;

  return true;
}

__attribute__ ((__deprecated__))
static bool bmm_dem_script_pushsediment(struct bmm_dem_opts* const opts,
    double const tspan, double const kcohes) {
  if (opts->script.n >= nmembof(opts->script.mode))
    return false;

  opts->script.tspan[opts->script.n] = tspan;
  opts->script.dt[opts->script.n] = dtstuff;
  opts->script.mode[opts->script.n] = BMM_DEM_MODE_SEDIMENT;
  opts->script.params[opts->script.n].sediment.kcohes = kcohes;

  ++opts->script.n;

  return true;
}

__attribute__ ((__nonnull__ (1, 2)))
static bool f(char const* const key, char const* const value,
    void* const ptr) {
  struct bmm_dem_opts* const opts = ptr;

  // TODO Refactor these at some point.

  if (strcmp(key, "script") == 0) {
    if (strcmp(value, "mix") == 0) {
      bmm_dem_script_pushidle(opts, 0.01);
      bmm_dem_script_pushcreate(opts, 0.75);
      bmm_dem_script_pushsediment(opts, 0.29, 1.0);
      bmm_dem_script_pushidle(opts, 0.20);
      // bmm_dem_script_pushidle(opts, 0.96);
    } else if (strcmp(value, "couple") == 0) {
    } else
      return false;
  } else if (strcmp(key, "verbose") == 0) {
    bool p;
    if (!bmm_str_strtob(&p, value))
      return false;

    opts->verbose = p;
  } else
    return false;

  return true;
}

__attribute__ ((__nonnull__))
int main(int const argc, char** const argv) {
  bmm_tle_reset(argv[0]);

#ifdef _GNU_SOURCE
#ifdef DEBUG
  int const excepts = feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  if (excepts == -1)
    BMM_TLE_STDS();
#endif
#endif

  struct bmm_dem_opts opts;
  bmm_dem_opts_def(&opts);

  if (!bmm_opt_parse((char const* const*) &argv[1], (size_t) (argc - 1),
        f, &opts)) {
    bmm_tle_put();

    return EXIT_FAILURE;
  }

  if (!bmm_dem_run_with(&opts)) {
    bmm_tle_put();

    return EXIT_FAILURE;
  }

#ifdef _GNU_SOURCE
#ifdef DEBUG
  if (feenableexcept(excepts) == -1)
    BMM_TLE_STDS();
#endif
#endif

  return EXIT_SUCCESS;
}
