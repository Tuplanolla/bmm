#ifndef BMM_GEOM_H
/// Geometry.
#define BMM_GEOM_H

#include <math.h>
#include <stddef.h>

#include "cpp.h"
#include "ext.h"
#include "fp.h"

/// The call `bmm_geom_ballvol(r, d)`
/// returns the volume of the `d`-dimensional ball of radius `r`.
/// The result is obtained by applying the recurrence relation
/// $V_0(r) = 1$, $V_1(r) = 2 r$, $V_d(r) = (\\twopi r^2 / d) V_{d - 2}(r)$.
__attribute__ ((__const__, __pure__))
inline double bmm_geom_ballvol(double const r, size_t d) {
  // return d == 0 ? 1.0 : d == 1 ? 2.0 * r :
  //   (M_2PI * BMM_POW(r, 2) / (double) d) * bmm_geom_ballvol(r, d - 2);

  double v = d % 2 == 0 ? 1.0 : 2.0 * r;

  while (d > 1) {
    v *= M_2PI * BMM_POW(r, 2) / (double) d;
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
  // return d == 0 ? 0.0 : d == 1 ? 2.0 : d == 2 ? M_2PI * r :
  //   (M_2PI * BMM_POW(r, 2) / (double) (d - 2)) * bmm_geom_ballsurf(r, d - 2);

  double a = d == 0 ? 0.0 : d % 2 == 1 ? 2.0 : M_2PI * r;

  while (d > 2) {
    d -= 2;
    a *= M_2PI * BMM_POW(r, 2) / (double) d;
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
  return bmm_geom_ballrmoi(d) * m * BMM_POW(r, 2);
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
  return bmm_geom_ballprmoi(d) * m * BMM_POW(r, 2);
}

/// The call `bmm_geom_ballmpd(d)`
/// returns the maximal packing density of balls in `d` dimensions.
/// The density is only known for
///
/// * 1,
/// * 2,
/// * 3,
/// * 8 and
/// * 24
///
/// dimensions; densities for all other dimensions are approximations.
__attribute__ ((__const__, __pure__))
inline double bmm_geom_ballmpd(size_t const d) {
  static double const delta[] = {
    NAN, 0.5, 0.28868, 0.17678, 0.125, 0.08839,
    0.07217, 0.0625, 0.0625, 0.04419, 0.03906, 0.03516,
    0.03704, 0.03516, 0.03608, 0.04419, 0.0625, 0.0625,
    0.07508, 0.08839, 0.13154, 0.17678, 0.33254, 0.5,
    1.0, 0.70711, 0.57735, 0.70711, 1.0, 0.70711,
    1.0, 1.2095, 2.5658, 2.2220, 2.2220, 2.8284, 4.4394
  };

  // TODO This.
  // return (pow(M_PI, (double) d / 2.0) / gsl_sf_gamma((double) d / 2.0 + 1)) * delta[d];

  switch (d) {
    case 1:
      return 1.0;
    case 2:
      return M_2PI / (4.0 * sqrt(2.0));
    case 3:
      return M_2PI / (6.0 * sqrt(2.0));
    case 8:
      return BMM_POW(M_2PI / 4.0, 4) / BMM_FACT(4);
    case 24:
      return BMM_POW(M_2PI / 2.0, 12) / BMM_FACT(12);
  }

  return NAN;
}

#endif
