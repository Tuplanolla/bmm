#ifndef BMM_GEOM_H
/// High-dimensional geometry.
#define BMM_GEOM_H

#include <math.h>
#include <stddef.h>

#include "cpp.h"
#include "ext.h"
#include "fp.h"
#include "size.h"

/// The call `bmm_geom_ballvol(r, d)`
/// returns the volume of the `d`-dimensional ball of radius `r`.
/// The result is obtained by applying the recurrence relation
/// $V_0(r) = 1$, $V_1(r) = 2 r$, $V_d(r) = (\\twopi r^2 / d) V_{d - 2}(r)$.
__attribute__ ((__const__, __pure__))
inline double bmm_geom_ballvol(double const r, size_t d) {
  // This would not be guaranteed to be tail-call optimized.
  // return d == 0 ? 1.0 : d == 1 ? 2.0 * r :
  //   (M_2PI * bmm_fp_pow(r, 2) / (double) d) * bmm_geom_ballvol(r, d - 2);

  double v = bmm_size_even(d) ? 1.0 : 2.0 * r;

  while (d > 1) {
    v *= M_2PI * bmm_fp_pow(r, 2) / (double) d;
    d -= 2;
  }

  return v;
}

/// The call `bmm_geom_ballsurf(r, d)`
/// returns the surface area of the `d`-dimensional ball of radius `r`.
/// The result is obtained by applying the recurrence relation
/// $A_0(r) = 0$, $A_1(r) = 2$, $A_2(r) = \\twopi r$,
/// $A_d(r) = (\\twopi r^2 / (d - 2)) A_{d - 2}(r)$.
__attribute__ ((__const__, __pure__))
inline double bmm_geom_ballsurf(double const r, size_t d) {
  // This would not be guaranteed to be tail-call optimized.
  // return d == 0 ? 0.0 : d == 1 ? 2.0 : d == 2 ? M_2PI * r :
  //   (M_2PI * bmm_fp_pow(r, 2) / (double) (d - 2)) * bmm_geom_ballsurf(r, d - 2);

  double a = d == 0 ? 0.0 : bmm_size_odd(d) ? 2.0 : M_2PI * r;

  while (d > 2) {
    d -= 2;
    a *= M_2PI * bmm_fp_pow(r, 2) / (double) d;
  }

  return a;
}

/// The call `bmm_geom_ballrmoi(d)`
/// returns the moment of inertia
/// of a `d`-dimensional homogeneous ball
/// as it rotates around any of its axes.
/// The result is obtained by applying the equation
/// $J_d(r) = ((d - 1) / (d + 2)) m r^2$.
__attribute__ ((__const__, __pure__))
inline double bmm_geom_ballrmoi(size_t const d) {
  if (d == 0)
    return (double) NAN;

  double const x = (double) d;

  return (x - 1.0) / (x + 2.0);
}

/// The call `bmm_geom_ballmoi(d)`
/// returns the moment of inertia
/// of a `d`-dimensional homogeneous ball with mass `m` and radius `r`
/// as it rotates around any of its axes.
/// See `bmm_geom_ballrmoi`.
__attribute__ ((__const__, __pure__))
inline double bmm_geom_ballmoi(double const m, double const r,
    size_t const d) {
  return bmm_geom_ballrmoi(d) * m * bmm_fp_pow(r, 2);
}

/// The call `bmm_geom_ballprmoi(d)`
/// returns the reduced moment of inertia
/// of a `d`-dimensional homogeneous ball
/// as it rotates around an axis perpendicular to all of its own axes.
/// This implies that the `d`-dimensional ball is embedded
/// into a `d + 1 + n` -dimensional space for some nonnegative `n`.
/// The result is obtained by applying the equation
/// $J_d(r) = (d / (d + 2)) m r^2$.
__attribute__ ((__const__, __pure__))
inline double bmm_geom_ballprmoi(size_t const d) {
  double const x = (double) d;

  return x / (x + 2.0);
}

/// The call `bmm_geom_ballpmoi(d)`
/// returns the moment of inertia
/// of a `d`-dimensional homogeneous ball with mass `m` and radius `r`
/// as it rotates around an axis perpendicular to all of its own axes.
/// See `bmm_geom_ballprmoi`.
__attribute__ ((__const__, __pure__))
inline double bmm_geom_ballpmoi(double const m, double const r,
    size_t const d) {
  return bmm_geom_ballprmoi(d) * m * bmm_fp_pow(r, 2);
}

/// The call `bmm_geom_ballmpcd(d)`
/// returns the maximal packing center density of balls in `d` dimensions.
/// The density is only known for
///
/// * the $1$-dimensional point lattice,
/// * the $2$-dimensional hexagonal lattice,
/// * the $3$-dimensional face-centered cubic lattice or
///   the equally dense hexagonal close-packing,
/// * the $8$-dimensional $E_8$ lattice and
/// * the $24$-dimensional Leech lattice.
///
/// Hypothetical answers exist for
///
/// * the $4$-dimensional $16$-cell honeycomb lattice,
/// * the $5$-dimensional $5$-demicubic honeycomb lattice,
/// * the $6$-dimensional $2_{22}$ honeycomb lattice and
/// * the $7$-dimensional $3_{31}$ honeycomb lattice.
///
/// Various approximate answers are used for other dimensions up to $36$.
/// Anything beyond $36$ dimensions is undefined.
/// The data source is `arXiv:math/0110009`.
__attribute__ ((__const__, __pure__))
inline double bmm_geom_ballmpcd(size_t const d) {
  static double const data[] = {
    1.0, 0.5, 0.28868, 0.17678, 0.125, 0.08839, 0.07217, 0.0625,
    0.0625, 0.04419, 0.03906, 0.03516, 0.03704, 0.03516, 0.03608, 0.04419,
    0.0625, 0.0625, 0.07508, 0.08839, 0.13154, 0.17678, 0.33254, 0.5,
    1.0, 0.70711, 0.57735, 0.70711, 1.0, 0.70711, 1.0, 1.2095,
    2.5658, 2.2220, 2.2220, 2.8284, 4.4394
  };

  return d < nmembof(data) ? data[d] : (double) NAN;
}

/// The call `bmm_geom_ballmpd(d)`
/// returns the maximal packing density of balls in `d` dimensions.
/// The density $\\Delta$ is related to the center density $\\delta$ by
/// $$\\frac \\Delta \\delta =
/// \\Big(\\frac \\twopi 2\\Big)^{d / 2} \\frac 1{\\Gamma(d / 2 + 1)}$$
/// or, in simplified form,
/// $$\\frac \\Delta \\delta =
/// \\Big(\\frac \\twopi 2\\Big)^{d / 2} \\frac 1{(d / 2)!}$$
/// for even $d$ and
/// $$\\frac \\Delta \\delta =
/// \\frac{(2 \\twopi)^{(d + 1) / 2}}{\\twopi / 2}
/// \\frac{((d + 1) / 2)!}{(d + 1)!} =
/// \\frac{\\twopi^{(d + 1) / 2}}{\\twopi / 2} \\frac 1{d!!}$$
/// for odd $d$.
///
/// See `bmm_geom_ballmpcd`.
__attribute__ ((__const__, __pure__))
inline double bmm_geom_ballmpd(size_t const d) {
  double c;

  if (bmm_size_even(d)) {
    size_t const n = d / 2;

    c = bmm_fp_pow(M_PI, n) / bmm_fp_fact(n, 1);
  } else {
    size_t const k = d + 1;
    size_t const n = k / 2;

    // This could be unstable.
    // c = (bmm_fp_pow(M_4PI, n) / M_PI) *
    //   (bmm_fp_fact(n, 1) / bmm_fp_fact(k, 1));
    c = (bmm_fp_pow(M_2PI, n) / M_PI) * (1.0 / bmm_fp_fact(d, 2));
  }

  return c * bmm_geom_ballmpcd(d);
}

#endif
