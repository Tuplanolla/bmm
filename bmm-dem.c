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

__attribute__ ((__nonnull__ (1, 2)))
static bool f(char const* const key, char const* const value,
    void* const ptr) {
  struct bmm_dem_opts* const opts = ptr;

  // TODO Refactor these at some point.
  double const dtstuff = 1.0e-4;
  size_t istage;

  opts->box.per[0] = true;
  opts->box.per[1] = false;

  // opts->cache.rcutoff /= 2.0;

  double const mu = 0.02;

  opts->part.rnew[0] = 2.0 * mu / (1.0 + sqrt(2.0));
  opts->part.rnew[1] = 4.0 * mu / (2.0 + sqrt(2.0));

  opts->part.rho = 1.0;
  opts->part.y = 1.0e+3;

  opts->comm.dt = 1.0e-3;

  if (strcmp(key, "script") == 0) {
    if (strcmp(value, "mix") == 0) {
      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_IDLE;
      opts->script.tspan[istage] = 0.005;
      opts->script.dt[istage] = dtstuff;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_CREATE;
      opts->script.params[istage].create.eta = 1.0;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_SEDIMENT;
      opts->script.tspan[istage] = 0.2;
      opts->script.dt[istage] = dtstuff;
      opts->script.params[istage].sediment.kcohes = 2.0;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_LINK;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_IDLE;
      opts->script.tspan[istage] = 0.2;
      opts->script.dt[istage] = dtstuff;
    } else if (strcmp(value, "couple") == 0) {
      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_TEST_COUPLE;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_IDLE;
      opts->script.tspan[istage] = 0.2;
      opts->script.dt[istage] = dtstuff;
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
