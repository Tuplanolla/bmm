#include <math.h>
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

  if (strcmp(key, "script") == 0) {
    // TODO Refactor these at some point.
    opts->box.x[1] = 0.1;
    opts->box.x[0] = opts->box.x[1];
    opts->box.per[0] = true;
    opts->box.per[1] = false;

    opts->part.rho = 2.7e+3;
    opts->part.ytens = 20.0e+9;
    opts->part.ycomp = 50.0e+9;
    opts->part.nu = 0.22;

    opts->part.strnew[0] = 1.0 - opts->gross.ds;
    opts->part.strnew[1] = 1.0 + opts->gross.ds;

    double const beta = 0.6; // From theory and experiments.
    double const sigmacrit = 280.0e+6; // From experiments.
    double const taucrit = beta * sigmacrit; // From theory.
    double const sigmacritt = sigmacrit;
    double const taucritt = beta * sigmacritt; // From theory.
    enum bmm_dem_tang const fric = opts->gross.fric;
    double ravg = 1.0e-3;
    double const eta = 1.0e+1; // Tar.
    double const eta2 = 1.0e+0; // Oil.
    double const eta3 = 1.0e-1; // Some other disgusting substance.
    double const a = 8.0e-9; // From theory and assumptions.
    double const mu = 0.7; // From experiments.
    double const kt = 1.0e+8; // From rolling and piling calibration tests.
    double const gammat = 1.0e+2; // From rolling and piling calibration tests.
    double const kn = 1.0e+7; // From theory and stress--strain tests.
    double const gamman = 1.0e+1; // From beam tests.
    double const barkn = 1.0e+7; // Equal to `kn`.
    double const bargamman = 1.0e+1; // Equal to `gamman`.
    double const barkt = 1.0e+7; // Small wrt `barkn`.
    double const bargammat = 1.0e+1; // Equal to `bargamman`.
    // What ought to hold.
    double const bshpp = (2.0 / 3.0) * (opts->part.ycomp /
        (1.0 - $(bmm_power, double)(opts->part.nu, 2)));
    double const more = bshpp * sqrt(ravg);
    double const knp = pow(M_PI * ravg * ravg *
        sigmacrit * more * more, 1.0 / 3.0);
    double const ap = gamman / more;
    // fprintf(stderr, "k^n > %g, A < %g\n", knp, ap);

    double const h = opts->gross.hfac;
    double const pdriv = sigmacrit * opts->gross.pfac;
    double const vdriv = 2.0;
    double const fadjust = 1.0e+8;
    double const padjust = 1.0e+8;

    double const rnew[] = {
      2.0 * ravg / (1.0 + sqrt(2.0)),
      4.0 * ravg / (2.0 + sqrt(2.0))
    };

    double dtstuff = 8.0e-7;
    size_t istage;

    if (strcmp(value, "fin") == 0) {
      dtstuff = 1.0e-7;
      opts->comm.dt = 1.0e-4;

      bmm_dem_opts_set_rnew(opts, rnew);

      opts->fault.njag = 4;
      opts->fault.hjag = h * 2.0 * bmm_ival_midpoint(opts->part.rnew);

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_CREATE_HC;
      opts->script.params[istage].create.eta = bmm_geom_ballmpd(BMM_NDIM);

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_PRESET0;
      opts->script.params[istage].preset.eta = eta;
      opts->script.params[istage].preset.kn = kn;
      opts->script.params[istage].preset.gamman = gamman;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_SEDIMENT;
      opts->script.tspan[istage] = 1.0e-3;
      opts->script.dt[istage] = dtstuff;
      opts->script.params[istage].sediment.kcohes = 2.0e+3;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_SEDIMENT;
      opts->script.tspan[istage] = 1.0e-3;
      opts->script.dt[istage] = dtstuff;
      opts->script.params[istage].sediment.kcohes = 0.2e+3;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_CLIP;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_EXPORT;
      opts->script.params[istage].expr.entropic = true;
      opts->script.params[istage].expr.str = "sediment";

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_LINK;

      dtstuff = 1.0e-8;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_PRESET1;
      opts->script.params[istage].preset.fric = fric;
      opts->script.params[istage].preset.eta2 = eta2;
      opts->script.params[istage].preset.a = a;
      opts->script.params[istage].preset.mu = mu;
      opts->script.params[istage].preset.kt = kt;
      opts->script.params[istage].preset.gammat = gammat;
      opts->script.params[istage].preset.kn = kn;
      opts->script.params[istage].preset.gamman = gamman;
      opts->script.params[istage].preset.barkn = barkn;
      opts->script.params[istage].preset.bargamman = bargamman;
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
      opts->script.mode[istage] = BMM_DEM_MODE_EXPORT;
      opts->script.params[istage].expr.entropic = false;
      opts->script.params[istage].expr.str = "fault";

      dtstuff = 1.0e-8;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_PRECRUNCH;
      opts->script.params[istage].precrunch.nlayer = 4.0;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_CRUNCH;
      opts->script.tspan[istage] = 0.5e-3;
      opts->script.dt[istage] = dtstuff;
      opts->script.params[istage].crunch.measure = false;
      opts->script.params[istage].crunch.v = 0.0;
      opts->script.params[istage].crunch.fadjust[0] = fadjust;
      opts->script.params[istage].crunch.fadjust[1] = padjust;
      opts->script.params[istage].crunch.p = pdriv;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_EXPORT;
      opts->script.params[istage].expr.entropic = false;
      opts->script.params[istage].expr.str = "shear";

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_ZEROEST;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_PRESET2;
      opts->script.params[istage].preset.eta3 = eta3;
      opts->script.params[istage].preset.sigmacrit = sigmacrit;
      opts->script.params[istage].preset.taucrit = taucrit;
      opts->script.params[istage].preset.sigmacritt = sigmacritt;
      opts->script.params[istage].preset.taucritt = taucritt;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_CRUNCH;
      opts->script.tspan[istage] = 50.0e-3;
      // opts->script.tspan[istage] = 10.0e-3; // This is good for testing.
      opts->script.dt[istage] = dtstuff;
      opts->script.params[istage].crunch.measure = true;
      opts->script.params[istage].crunch.v = vdriv;
      opts->script.params[istage].crunch.fadjust[0] = fadjust;
      opts->script.params[istage].crunch.fadjust[1] = padjust;
      opts->script.params[istage].crunch.p = pdriv;

      // Strain 100 %.
      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_EXPORT;
      opts->script.params[istage].expr.entropic = false;
      opts->script.params[istage].expr.str = "fragment";
    } else if (strcmp(value, "shear") == 0) {
      opts->box.x[1] = 0.025;
      opts->box.x[0] = opts->box.x[1];

      bmm_dem_opts_set_rnew(opts, rnew);

      dtstuff = 2.0e-7;
      opts->comm.dt = 1.0e-4;

      opts->fault.njag = 2;
      opts->fault.hjag = h * 2.0 * bmm_ival_midpoint(opts->part.rnew);

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_CREATE_HC;
      opts->script.params[istage].create.eta = bmm_geom_ballmpd(BMM_NDIM);

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_PRESET0;
      opts->script.params[istage].preset.eta = eta;
      opts->script.params[istage].preset.kn = kn;
      opts->script.params[istage].preset.gamman = gamman;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_SEDIMENT;
      opts->script.tspan[istage] = 1.0e-3;
      opts->script.dt[istage] = dtstuff;
      opts->script.params[istage].sediment.kcohes = 4.0e+3;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_SEDIMENT;
      opts->script.tspan[istage] = 1.0e-3;
      opts->script.dt[istage] = dtstuff;
      opts->script.params[istage].sediment.kcohes = 0.4e+3;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_CLIP;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_EXPORT;
      opts->script.params[istage].expr.entropic = true;
      opts->script.params[istage].expr.str = "sediment";

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_LINK;

      dtstuff = 2.0e-8;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_PRESET1;
      opts->script.params[istage].preset.fric = fric;
      opts->script.params[istage].preset.eta2 = eta2;
      opts->script.params[istage].preset.a = a;
      opts->script.params[istage].preset.mu = mu;
      opts->script.params[istage].preset.kt = kt;
      opts->script.params[istage].preset.gammat = gammat;
      opts->script.params[istage].preset.kn = kn;
      opts->script.params[istage].preset.gamman = gamman;
      opts->script.params[istage].preset.barkn = barkn;
      opts->script.params[istage].preset.bargamman = bargamman;
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
      opts->script.mode[istage] = BMM_DEM_MODE_EXPORT;
      opts->script.params[istage].expr.entropic = false;
      opts->script.params[istage].expr.str = "fault";

      dtstuff = 2.0e-8;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_PRECRUNCH;
      opts->script.params[istage].precrunch.nlayer = 2.0;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_CRUNCH;
      opts->script.tspan[istage] = 0.5e-3;
      opts->script.dt[istage] = dtstuff;
      opts->script.params[istage].crunch.measure = false;
      opts->script.params[istage].crunch.v = 0.0;
      opts->script.params[istage].crunch.fadjust[0] = fadjust;
      opts->script.params[istage].crunch.fadjust[1] = padjust;
      opts->script.params[istage].crunch.p = pdriv;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_EXPORT;
      opts->script.params[istage].expr.entropic = false;
      opts->script.params[istage].expr.str = "shear";

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_ZEROEST;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_PRESET2;
      opts->script.params[istage].preset.eta3 = eta3;
      opts->script.params[istage].preset.sigmacrit = sigmacrit;
      opts->script.params[istage].preset.taucrit = taucrit;
      opts->script.params[istage].preset.sigmacritt = sigmacritt;
      opts->script.params[istage].preset.taucritt = taucritt;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_CRUNCH;
      opts->script.tspan[istage] = 7.0e-3;
      opts->script.dt[istage] = dtstuff;
      opts->script.params[istage].crunch.measure = true;
      opts->script.params[istage].crunch.v = vdriv;
      opts->script.params[istage].crunch.fadjust[0] = fadjust;
      opts->script.params[istage].crunch.fadjust[1] = padjust;
      opts->script.params[istage].crunch.p = pdriv;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_EXPORT;
      opts->script.params[istage].expr.entropic = false;
      opts->script.params[istage].expr.str = "fragment";
    } else if (strcmp(value, "beam") == 0) {
      dtstuff = 4.0e-8;
      opts->comm.dt = 4.0e-5;

      double const rnew[] = {2.0e-3, 2.0e-3 + 1.0e-6};
      bmm_dem_opts_set_rnew(opts, rnew);

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_IDLE;
      opts->script.tspan[istage] = 0.1e-3;
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
      opts->script.params[istage].preset.eta = eta;
      opts->script.params[istage].preset.kn = kn;
      opts->script.params[istage].preset.gamman = gamman;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_PRESET1;
      opts->script.params[istage].preset.fric = fric;
      opts->script.params[istage].preset.eta2 = eta2;
      opts->script.params[istage].preset.a = a;
      opts->script.params[istage].preset.mu = mu;
      opts->script.params[istage].preset.kt = kt;
      opts->script.params[istage].preset.gammat = gammat;
      opts->script.params[istage].preset.kn = kn;
      opts->script.params[istage].preset.gamman = gamman;
      opts->script.params[istage].preset.barkn = barkn;
      opts->script.params[istage].preset.bargamman = bargamman;
      opts->script.params[istage].preset.barkt = barkt;
      opts->script.params[istage].preset.bargammat = bargammat;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_PRESET2;
      opts->script.params[istage].preset.eta3 = eta3;
      opts->script.params[istage].preset.sigmacrit = sigmacrit;
      opts->script.params[istage].preset.taucrit = taucrit;
      opts->script.params[istage].preset.sigmacritt = sigmacritt;
      opts->script.params[istage].preset.taucritt = taucritt;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_GRAVY;
      opts->script.tspan[istage] = 10.0e-3;
      opts->script.dt[istage] = dtstuff * 0.5;
      opts->script.params[istage].gravy.g = 3.0e+5;
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
      opts->comm.dt = 4.0e-5;
      dtstuff = 1.0e-7;

      double const rnew[] = {
        3.0e-3,
        3.0e-3 * 1.07
      };
      bmm_dem_opts_set_rnew(opts, rnew);

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_PRESET0;
      opts->script.params[istage].preset.eta = eta;
      opts->script.params[istage].preset.kn = kn;
      opts->script.params[istage].preset.gamman = gamman;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_PRESET1;
      opts->script.params[istage].preset.fric = fric;
      opts->script.params[istage].preset.eta2 = eta2;
      opts->script.params[istage].preset.a = a;
      opts->script.params[istage].preset.mu = mu;
      opts->script.params[istage].preset.kt = kt;
      opts->script.params[istage].preset.gammat = gammat;
      opts->script.params[istage].preset.kn = kn;
      opts->script.params[istage].preset.gamman = gamman;
      opts->script.params[istage].preset.barkn = barkn;
      opts->script.params[istage].preset.bargamman = bargamman;
      opts->script.params[istage].preset.barkt = barkt;
      opts->script.params[istage].preset.bargammat = bargammat;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_PRESET2;
      opts->script.params[istage].preset.eta3 = eta3;
      opts->script.params[istage].preset.sigmacrit = sigmacrit;
      opts->script.params[istage].preset.taucrit = taucrit;
      opts->script.params[istage].preset.sigmacritt = sigmacritt;
      opts->script.params[istage].preset.taucritt = taucritt;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_CREATE_PLANE;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_CREATE_COUPLE;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_GRAVY;
      opts->script.tspan[istage] = 20.0e-3;
      opts->script.dt[istage] = dtstuff;
      opts->script.params[istage].gravy.g = 3.0e+4;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_LINK;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_IDLE;
      opts->script.tspan[istage] = 20.0e-3;
      opts->script.dt[istage] = dtstuff;
    } else if (strcmp(value, "pile") == 0) {
      opts->comm.dt = 4.0e-5;
      dtstuff = 1.0e-8;

      double const rnew[] = {
        2.0e-3,
        2.0e-3 + 2.0e-5
      };
      bmm_dem_opts_set_rnew(opts, rnew);

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_IDLE;
      opts->script.tspan[istage] = 0.1e-3;
      opts->script.dt[istage] = dtstuff;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_PRESET0;
      opts->script.params[istage].preset.eta = eta;
      opts->script.params[istage].preset.kn = kn;
      opts->script.params[istage].preset.gamman = gamman;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_PRESET1;
      opts->script.params[istage].preset.fric = fric;
      opts->script.params[istage].preset.eta2 = eta2;
      opts->script.params[istage].preset.a = a;
      opts->script.params[istage].preset.mu = mu;
      opts->script.params[istage].preset.kt = kt;
      opts->script.params[istage].preset.gammat = gammat;
      opts->script.params[istage].preset.kn = kn;
      opts->script.params[istage].preset.gamman = gamman;
      opts->script.params[istage].preset.barkn = barkn;
      opts->script.params[istage].preset.bargamman = bargamman;
      opts->script.params[istage].preset.barkt = barkt;
      opts->script.params[istage].preset.bargammat = bargammat;

      istage = bmm_dem_script_addstage(opts);
      opts->script.mode[istage] = BMM_DEM_MODE_PRESET2;
      opts->script.params[istage].preset.eta3 = eta3;
      opts->script.params[istage].preset.sigmacrit = sigmacrit;
      opts->script.params[istage].preset.taucrit = taucrit;
      opts->script.params[istage].preset.sigmacritt = sigmacritt;
      opts->script.params[istage].preset.taucritt = taucritt;

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
      opts->script.params[istage].gravy.g = 3.0e+4;

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
  } else if (strcmp(key, "fric") == 0) {
    if (strcmp(value, "cs") == 0)
      opts->gross.fric = BMM_DEM_TANG_CS;
    else if (strcmp(value, "hw") == 0)
      opts->gross.fric = BMM_DEM_TANG_HW;
    else if (strcmp(value, "none") == 0)
      opts->gross.fric = BMM_DEM_TANG_NONE;
    else
      return false;
  } else if (strcmp(key, "hfac") == 0) {
    double x;
    if (!bmm_str_strtod(&x, value))
      return false;

    opts->gross.hfac = x;
  } else if (strcmp(key, "pfac") == 0) {
    double x;
    if (!bmm_str_strtod(&x, value))
      return false;

    opts->gross.pfac = x;
  } else if (strcmp(key, "ds") == 0) {
    double x;
    if (!bmm_str_strtod(&x, value))
      return false;

    opts->gross.ds = x;
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
