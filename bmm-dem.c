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

__attribute__ ((__nonnull__))
static bool ind_straight(double const *const x, void *const cls) {
  struct bmm_dem_opts *const opts = cls;

  return x[1] > opts->box.x[1] / 2.0;
}

__attribute__ ((__nonnull__))
static bool ind_wave(double const *const x, void *const cls) {
  struct bmm_dem_opts *const opts = cls;

  return x[1] > opts->box.x[1] / 2.0 + (opts->fault.hjag / 2.0) *
    sin(M_2PI * (double) opts->fault.njag * x[0] / opts->box.x[0]);
}

__attribute__ ((__nonnull__ (1, 2)))
static bool f(char const *const key, char const *const value,
    void *const ptr) {
  struct bmm_dem_opts *const opts = ptr;

  // TODO Refactor these at some point.
  opts->box.x[1] = 0.1;
  opts->box.x[0] = opts->box.x[1];
  opts->box.per[0] = true;
  opts->box.per[1] = false;

  double s = 0.3;
  opts->part.strnew[0] = 1.0 - s;
  opts->part.strnew[1] = 1.0 + s;

  opts->part.rho = 2.65e+3;
  opts->part.ytens = 17.0e+9;
  opts->part.ycomp = 51.0e+9;
  opts->part.nu = 0.215;

  double const beta = 0.5;
  double const sigmacrit = 283.0e+6; // From experiments.
  double const taucrit = beta * sigmacrit; // From theory.
  // double const sigmacritt = 10.0e+6; // From experiments.
  double const sigmacritt = sigmacrit; // From experiments.
  double const taucritt = beta * sigmacritt; // From theory.
  bool const statfric = true;
  double ravg = 1.0e-3;
  double const eta = 1.0e+3; // Water.
  double const eta2 = 2.0e+3; // Glue.
  double const a = 7.6e-3; // From theory.
  double const mu = 0.725; // From experiments.
  double const kt = 1.0e+6;
  double const gammat = 1.0e+1;
  double const kn = 3.0e+7;
  double const gamman = 1.0e+2; // Critical damping.
  double const barkt = kn * 0.1; // Beam has neutral axis in the middle.
  double const bargammat = 1.0e+2; // Critical damping.

  double const pdriv = 183.0e+6;
  double const vdriv = 4.0;
  double const fadjust = 4.0e+9;

  double dtstuff = 8.0e-7;
  size_t istage;

  if (strcmp(key, "script") == 0) {
    if (strcmp(value, "fin") == 0) {
      dtstuff = 1.0e-7;
      opts->comm.dt = 1.0e-4;

      double const rnew[] = {
        2.0 * ravg / (1.0 + sqrt(2.0)),
        4.0 * ravg / (2.0 + sqrt(2.0))
      };
      bmm_dem_opts_set_rnew(opts, rnew);

      opts->fault.njag = 4;
      opts->fault.hjag = 4.0 * 2.0 * bmm_ival_midpoint(opts->part.rnew);

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_IDLE;
      opts->script.tspan[istage] = 0.03e-3;
      opts->script.dt[istage] = dtstuff;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_CREATE_HC;
      opts->script.params[istage].create.eta = bmm_geom_ballmpd(BMM_NDIM);

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_PRESET0;
      opts->script.params[istage].preset.statfric = statfric;
      opts->script.params[istage].preset.eta = eta;
      opts->script.params[istage].preset.kn = kn;
      opts->script.params[istage].preset.gamman = gamman;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_SEDIMENT;
      opts->script.tspan[istage] = 2.0e-3;
      opts->script.dt[istage] = dtstuff;
      opts->script.params[istage].sediment.kcohes = 0.5e+4;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_SEDIMENT;
      opts->script.tspan[istage] = 1.0e-3;
      opts->script.dt[istage] = dtstuff;
      opts->script.params[istage].sediment.kcohes = 1.0e+3;

      dtstuff = 1.0e-8;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_CLIP;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_LINK;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_PRESET1;
      opts->script.params[istage].preset.eta2 = eta2;
      opts->script.params[istage].preset.a = a;
      opts->script.params[istage].preset.mu = mu;
      opts->script.params[istage].preset.kt = kt;
      opts->script.params[istage].preset.gammat = gammat;
      opts->script.params[istage].preset.kn = kn;
      opts->script.params[istage].preset.gamman = gamman;
      opts->script.params[istage].preset.barkt = barkt;
      opts->script.params[istage].preset.bargammat = bargammat;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_SET_DENSITY;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_FAULT;
      opts->script.params[istage].fault.ind = ind_wave;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_IDLE;
      opts->script.tspan[istage] = 0.5e-3;
      opts->script.dt[istage] = dtstuff;

      dtstuff = 10.0e-9;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_PRESET2;
      opts->script.params[istage].preset.sigmacrit = sigmacrit;
      opts->script.params[istage].preset.taucrit = taucrit;
      opts->script.params[istage].preset.sigmacritt = sigmacritt;
      opts->script.params[istage].preset.taucritt = taucritt;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_CRUNCH;
      opts->script.tspan[istage] = 8.0e-3;
      // opts->script.tspan[istage] = 20.0e-3;
      opts->script.dt[istage] = dtstuff;
      opts->script.params[istage].crunch.nlayer = 3.0;
      opts->script.params[istage].crunch.v = vdriv;
      opts->script.params[istage].crunch.fadjust[0] =
        opts->script.params[istage].crunch.fadjust[1] = fadjust;
      opts->script.params[istage].crunch.p = pdriv;
    } else if (strcmp(value, "shear") == 0) {
      double ravg = 4.0e-3;
      dtstuff = 6.0e-7;
      // dtstuff = 4.0e-9;
      opts->comm.dt = 1.0e-4;

      double const rnew[] = {
        2.0 * ravg / (1.0 + sqrt(2.0)),
        4.0 * ravg / (2.0 + sqrt(2.0))
      };
      bmm_dem_opts_set_rnew(opts, rnew);

      opts->fault.njag = 2;
      opts->fault.hjag = 2.0 * 2.0 * bmm_ival_midpoint(opts->part.rnew);

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_IDLE;
      opts->script.tspan[istage] = 0.03e-3;
      opts->script.dt[istage] = dtstuff;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_CREATE_HC;
      opts->script.params[istage].create.eta = bmm_geom_ballmpd(BMM_NDIM);

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_PRESET0;
      opts->script.params[istage].preset.statfric = statfric;
      opts->script.params[istage].preset.eta = eta;
      opts->script.params[istage].preset.kn = kn;
      opts->script.params[istage].preset.gamman = gamman;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_SEDIMENT;
      opts->script.tspan[istage] = 4.5e-3;
      opts->script.dt[istage] = dtstuff;
      opts->script.params[istage].sediment.kcohes = 3.0e+4;

      dtstuff = 6.0e-8;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_SEDIMENT;
      opts->script.tspan[istage] = 1.0e-3;
      opts->script.dt[istage] = dtstuff;
      opts->script.params[istage].sediment.kcohes = 1.0e+4;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_CLIP;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_LINK;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_PRESET1;
      opts->script.params[istage].preset.eta2 = eta2;
      opts->script.params[istage].preset.a = a;
      opts->script.params[istage].preset.mu = mu;
      opts->script.params[istage].preset.kt = kt;
      opts->script.params[istage].preset.gammat = gammat;
      opts->script.params[istage].preset.kn = kn;
      opts->script.params[istage].preset.gamman = gamman;
      opts->script.params[istage].preset.barkt = barkt;
      opts->script.params[istage].preset.bargammat = bargammat;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_SET_DENSITY;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_FAULT;
      opts->script.params[istage].fault.ind = ind_wave;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_IDLE;
      opts->script.tspan[istage] = 0.5e-3;
      opts->script.dt[istage] = dtstuff;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_PRESET2;
      opts->script.params[istage].preset.sigmacrit = sigmacrit;
      opts->script.params[istage].preset.taucrit = taucrit;
      opts->script.params[istage].preset.sigmacritt = sigmacritt;
      opts->script.params[istage].preset.taucritt = taucritt;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_CRUNCH;
      opts->script.tspan[istage] = 8.0e-3;
      // opts->script.tspan[istage] = 20.0e-3;
      opts->script.dt[istage] = dtstuff;
      opts->script.params[istage].crunch.nlayer = 2.0;
      opts->script.params[istage].crunch.v = vdriv;
      opts->script.params[istage].crunch.fadjust[0] =
        opts->script.params[istage].crunch.fadjust[1] = fadjust;
      opts->script.params[istage].crunch.p = pdriv;
    } else if (strcmp(value, "beam") == 0) {
      dtstuff = 2.0e-7;
      opts->comm.dt = 4.0e-5;

      double const rnew[] = {2.0e-3, 2.0e-3 + 1.0e-6};
      bmm_dem_opts_set_rnew(opts, rnew);

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_IDLE;
      opts->script.tspan[istage] = 0.03e-3;
      opts->script.dt[istage] = dtstuff;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_CREATE_BEAM;
      opts->script.params[istage].test.eta = bmm_geom_ballmpd(BMM_NDIM);
      opts->script.params[istage].test.layers = 24.0;
      opts->script.params[istage].test.slices = 42.0;
      // opts->script.params[istage].test.layers = 16.0;
      // opts->script.params[istage].test.slices = 16.0;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_LINK;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_PRESET0;
      opts->script.params[istage].preset.statfric = statfric;
      opts->script.params[istage].preset.eta = eta;
      opts->script.params[istage].preset.kn = kn;
      opts->script.params[istage].preset.gamman = gamman;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_PRESET1;
      opts->script.params[istage].preset.eta2 = eta2;
      opts->script.params[istage].preset.a = a;
      opts->script.params[istage].preset.mu = mu;
      opts->script.params[istage].preset.kt = kt;
      opts->script.params[istage].preset.gammat = gammat;
      opts->script.params[istage].preset.kn = kn;
      opts->script.params[istage].preset.gamman = gamman;
      opts->script.params[istage].preset.barkt = barkt;
      opts->script.params[istage].preset.bargammat = bargammat;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_PRESET2;
      opts->script.params[istage].preset.sigmacrit = sigmacrit;
      opts->script.params[istage].preset.taucrit = taucrit;
      opts->script.params[istage].preset.sigmacritt = sigmacritt;
      opts->script.params[istage].preset.taucritt = taucritt;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_GRAVY;
      opts->script.tspan[istage] = 10.0e-3;
      opts->script.dt[istage] = dtstuff * 0.5;
      opts->script.params[istage].gravy.g = 3.0e+4;
    } else if (strcmp(value, "mix") == 0) {
      double const mu = 2.0e-3;

      double const rnew[] = {
        2.0 * mu / (1.0 + sqrt(2.0)),
        4.0 * mu / (2.0 + sqrt(2.0))
      };
      bmm_dem_opts_set_rnew(opts, rnew);

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_IDLE;
      opts->script.tspan[istage] = 0.03e-3;
      opts->script.dt[istage] = dtstuff;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_CREATE_HC;
      opts->script.params[istage].create.eta = bmm_geom_ballmpd(BMM_NDIM);

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_SEDIMENT;
      opts->script.tspan[istage] = 1.5e-3;
      opts->script.dt[istage] = dtstuff;
      opts->script.params[istage].sediment.kcohes = 3.0e+4;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_CLIP;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_LINK;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_IDLE;
      opts->script.tspan[istage] = 6.0e-3;
      opts->script.dt[istage] = dtstuff;
    } else if (strcmp(value, "roll") == 0) {
      dtstuff = 2.0e-7;

      double const rnew[] = {
        2.078e-3,
        2.078e-3 * 1.07
      };
      bmm_dem_opts_set_rnew(opts, rnew);

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_IDLE;
      opts->script.tspan[istage] = 0.03e-3;
      opts->script.dt[istage] = dtstuff;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_CREATE_PLANE;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_CREATE_COUPLE;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_GRAVY;
      opts->script.tspan[istage] = 3.0e-3;
      opts->script.dt[istage] = dtstuff;
      opts->script.params[istage].gravy.g = 3.0e+4;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_LINK;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_IDLE;
      opts->script.tspan[istage] = 6.0e-3;
      opts->script.dt[istage] = dtstuff;
    } else if (strcmp(value, "pile") == 0) {
      dtstuff = 2.0e-8;

      double const rnew[] = {
        2.078e-3,
        2.078e-3 * 1.37
      };
      bmm_dem_opts_set_rnew(opts, rnew);

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_IDLE;
      opts->script.tspan[istage] = 0.03e-3;
      opts->script.dt[istage] = dtstuff;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_CREATE_PILE;
      opts->script.params[istage].test.eta = bmm_geom_ballmpd(BMM_NDIM);
      /*
      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_CREATE_HC;
      opts->script.params[istage].create.eta = 0.5;
      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_CREATE_PLANE;
      */

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_GRAVY;
      opts->script.tspan[istage] = 30.0e-3;
      opts->script.dt[istage] = dtstuff;
      opts->script.params[istage].gravy.g = 3.0e+5;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_LINK;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_IDLE;
      opts->script.tspan[istage] = 3.0e-3;
      opts->script.dt[istage] = dtstuff;
    } else if (strcmp(value, "couple") == 0) {
      dtstuff = 2.0e-8;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_IDLE;
      opts->script.tspan[istage] = 0.06e-3;
      opts->script.dt[istage] = dtstuff;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_CREATE_COUPLE;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_IDLE;
      opts->script.tspan[istage] = 0.6e-3;
      opts->script.dt[istage] = dtstuff;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_LINK;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_IDLE;
      opts->script.tspan[istage] = 6.0e-3;
      opts->script.dt[istage] = dtstuff;
    } else if (strcmp(value, "triplet") == 0) {
      dtstuff = 2.0e-8;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_IDLE;
      opts->script.tspan[istage] = 0.06e-3;
      opts->script.dt[istage] = dtstuff;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_CREATE_TRIPLET;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_IDLE;
      opts->script.tspan[istage] = 0.6e-3;
      opts->script.dt[istage] = dtstuff;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_LINK;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_IDLE;
      opts->script.tspan[istage] = 6.0e-3;
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
