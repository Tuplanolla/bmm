#ifndef BMM_GEOM_H
/// Geometry.
#define BMM_GEOM_H

#include <math.h>
#include <stddef.h>

#include "ext.h"
#include "fp.h"

/// The call `bmm_geom_ballvol(r, d)`
/// returns the volume of the `d`-dimensional ball of radius `r`.
/// The result is obtained by applying the recurrence relation
/// $V_0(r) = 1$, $V_1(r) = 2 r$, $V_d(r) = (\\twopi r^2 / d) V_{d - 2}(r)$.
__attribute__ ((__const__, __pure__))
inline double bmm_geom_ballvol(double const r, size_t d) {
  // return d == 0 ? 1.0 : d == 1 ? 2.0 * r :
  //   (M_2PI * bmm_fp_sq(r) / (double) d) * bmm_geom_ballvol(r, d - 2);

  double v = d % 2 == 0 ? 1.0 : 2.0 * r;

  while (d > 1) {
    v *= M_2PI * bmm_fp_sq(r) / (double) d;
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
  //   (M_2PI * bmm_fp_sq(r) / (double) (d - 2)) * bmm_geom_ballsurf(r, d - 2);

  double a = d == 0 ? 0.0 : d % 2 == 1 ? 2.0 : M_2PI * r;

  while (d > 2) {
    d -= 2;
    a *= M_2PI * bmm_fp_sq(r) / (double) d;
  }

  return a;
}

/// The call `bmm_geom_ballmoi(r, d)`
/// returns the moment of inertia
/// of a `d`-dimensional homogeneous unit mass ball of radius `r`
/// as it rotates around any of its axes.
/// To obtain the moment of inertia of a ball with mass `m`,
/// multiply the return value with `m`.
/// The result is obtained by applying the equation
/// $J_d(r) = ((d - 1) / (d + 2)) m r^2$.
__attribute__ ((__const__, __pure__))
inline double bmm_geom_ballmoi(double const r, size_t const d) {
  if (d == 0)
    return (double) NAN;

  double const x = (double) d;

  return ((x - 1.0) / (x + 2.0)) * bmm_fp_sq(r);
}

/// The call `bmm_geom_ballpmoi(r, d)`
/// returns the moment of inertia
/// of a `d`-dimensional homogeneous unit mass ball of radius `r`
/// as it rotates around an axis perpendicular to all of its own axes.
/// This implies that the `d`-dimensional ball is embedded
/// into a `d + 1 + n` -dimensional space for some nonnegative `n`.
/// To obtain the moment of inertia of a ball with mass `m`,
/// multiply the return value with `m`.
/// The result is obtained by applying the equation
/// $J_d(r) = (d / (d + 2)) m r^2$.
__attribute__ ((__const__, __pure__))
inline double bmm_geom_ballpmoi(double const r, size_t const d) {
  double const x = (double) d;

  return (x / (x + 2.0)) * bmm_fp_sq(r);
}

#endif
