// Geometry.
#ifndef BMM_GEOM_H
#define BMM_GEOM_H

#include "ext.h"
#include "fp.h"
#include <math.h>
#include <stddef.h>

// The call `bmm_geom_ballvol(r, d)` returns the volume
// of the `d`-dimensional ball of radius `r`.
// The result is obtained by applying the recurrence relation
// $V_0(r) = 1$, $V_1(r) = 2 r$, $V_d(r) = (\twopi r^2 / d) V_{d - 2}(r)$.
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

// The call `bmm_geom_ballsurf(r, d)` returns the surface area
// of the `d`-dimensional ball of radius `r`.
// The result is obtained by applying the recurrence relation
// $A_0(r) = 0$, $A_1(r) = 2$, $A_2(r) = \twopi r$,
// $A_d(r) = (\twopi r^2 / (d - 2)) A_{d - 2}(r)$.
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

// The call `bmm_geom_ballmoi(r, d)` returns the moment of inertia
// of the `d`-dimensional solid unit mass ball of radius `r`.
// The result is obtained by evaluating the integral
// $J_d(r) = ...$.
__attribute__ ((__const__, __pure__))
inline double bmm_geom_ballmoi(double const r, size_t const d) {
  switch (d) {
    case 1:
      return 0.0;
    case 2:
      return 0.5 * bmm_fp_sq(r);
    case 3:
      return 0.4 * bmm_fp_sq(r);
    default:
      // TODO Use a fallback.
      return NAN;
  }
}

#endif
