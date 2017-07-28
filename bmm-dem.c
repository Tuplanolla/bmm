#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "dem.h"
#include "ext.h"
#include "fp.h"
#include "geom.h"
#include "opt.h"
#include "str.h"
#include "tle.h"

__attribute__ ((__nonnull__ (1, 2)))
static bool f(char const *const key, char const *const value,
    void *const ptr) {
  struct bmm_dem_opts *const opts = ptr;

  // TODO Refactor these at some point.
  double const dtstuff = 1.0e-6;
  size_t istage;

  opts->box.x[0] = 0.1;
  opts->box.x[1] = 0.1;
  opts->box.per[0] = true;
  opts->box.per[1] = false;

  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    opts->cache.ncell[idim] = 12;

  opts->cache.rcutoff = (double) INFINITY;
  for (size_t idim = 0; idim < BMM_NDIM; ++idim)
    opts->cache.rcutoff = fmin(opts->cache.rcutoff,
        opts->box.x[idim] / (double) (opts->cache.ncell[idim] - 2));

  double const mu = 2.0e-3;

  opts->part.rnew[0] = 2.0 * mu / (1.0 + sqrt(2.0));
  opts->part.rnew[1] = 4.0 * mu / (2.0 + sqrt(2.0));

  // Static pile of sand test.
  opts->part.rnew[0] = 2.078e-3;
  opts->part.rnew[1] = opts->part.rnew[0] + 1.0e-9;

  // Now in SI base units!
  opts->part.rho = 2.7e+3;
  // TODO For granite this should be closer to `e+9`,
  // but that explodes with this large `dt`.
  opts->part.y = 52.0e+8;
  opts->part.nu = 0.2;

  opts->comm.dt = 2.0e-5;

  if (strcmp(key, "script") == 0) {
    if (strcmp(value, "mix") == 0) {
      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_IDLE;
      opts->script.tspan[istage] = 0.03e-3;
      opts->script.dt[istage] = dtstuff;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_CREATE;
      opts->script.params[istage].create.eta = bmm_geom_ballmpd(BMM_NDIM);

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_SEDIMENT;
      opts->script.tspan[istage] = 3.0e-3;
      opts->script.dt[istage] = dtstuff;
      opts->script.params[istage].sediment.kcohes = 3.0e+4;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_CLIP;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_LINK;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_IDLE;
      opts->script.tspan[istage] = 3.0e-3;
      opts->script.dt[istage] = dtstuff;
    } else if (strcmp(value, "test") == 0) {
      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_CREATE;
      opts->script.params[istage].create.eta = bmm_geom_ballmpd(BMM_NDIM);

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_SEDIMENT;
      opts->script.tspan[istage] = 2.0e-3;
      opts->script.dt[istage] = dtstuff;
      opts->script.params[istage].sediment.kcohes = 3.0e+4;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_CLIP;
    } else if (strcmp(value, "couple") == 0) {
      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_IDLE;
      opts->script.tspan[istage] = 0.06e-3;
      opts->script.dt[istage] = dtstuff;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_TEST_COUPLE;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_IDLE;
      opts->script.tspan[istage] = 3.0e-3;
      opts->script.dt[istage] = dtstuff;
    } else
      return false;
  } else if (strcmp(key, "verbose") == 0) {
    bool p;
    if (!bmm_str_strtob(&p, value))
      return false;

    opts->verbose = p;
  } else if (strcmp(key, "trap") == 0) {
    bool p;
    if (!bmm_str_strtob(&p, value))
      return false;

    opts->trap.enabled = p;
  } else
    return false;

  return true;
}

__attribute__ ((__nonnull__))
int main(int const argc, char **const argv) {
  bmm_tle_reset(argv[0]);

  struct bmm_dem_opts opts;
  bmm_dem_opts_def(&opts);

  if (!bmm_opt_parse((char const *const *) &argv[1], (size_t) (argc - 1),
        f, &opts)) {
    bmm_tle_put();

    return EXIT_FAILURE;
  }

  if (!bmm_dem_run_with(&opts)) {
    bmm_tle_put();

    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
