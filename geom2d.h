// Planar geometry.
#ifndef BMM_GEOM2D_H
#define BMM_GEOM2D_H

#include "ext.h"
#include "fp.h"
#include <math.h>
#include <stddef.h>

// The call `bmm_geom2d_norm2(r)` returns the length $r^2$
// of the two-dimensional vector `r`.
__attribute__ ((__pure__))
inline double bmm_geom2d_norm2(double const* const r) {
  double d = 0.0;
  for (size_t idim = 0; idim < 2; ++idim)
    d += bmm_fp_sq(r[idim]);

  return d;
}

// The call `bmm_geom2d_norm(r)` returns the length $r$
// of the two-dimensional vector `r`.
__attribute__ ((__pure__))
inline double bmm_geom2d_norm(double const* const r) {
  return sqrt(bmm_geom2d_norm2(r));
}

// The call `bmm_geom2d_angle(r)` returns the signed angle $\phi$
// of the two-dimensional vector `r`.
__attribute__ ((__pure__))
inline double bmm_geom2d_angle(double const* const r) {
  return atan2(r[1], r[0]);
}

// The call `bmm_geom2d_to_polar2(polar2, cart)` maps
// the two-dimensional cartesian coordinates `cart` of the form $(x, y)$
// to the corresponding polar coordinates `polar2` of the form $(r^2, \phi)$.
inline void bmm_geom2d_to_polar2(double* const polar2,
    double const* const cart) {
  polar2[0] = bmm_geom2d_norm2(cart);
  polar2[1] = bmm_geom2d_angle(cart);
}

// The call `bmm_geom2d_to_polar(polar, cart)` maps
// the two-dimensional cartesian coordinates `cart` of the form $(x, y)$
// to the corresponding polar coordinates `polar` of the form $(r, \phi)$.
inline void bmm_geom2d_to_polar(double* const polar,
    double const* const cart) {
  polar[0] = bmm_geom2d_norm(cart);
  polar[1] = bmm_geom2d_angle(cart);
}

// The call `bmm_geom2d_from_polar(cart, polar)` maps
// the two-dimensional polar coordinates `polar` of the form $(r, \phi)$
// to the corresponding cartesian coordinates `cart` of the form $(x, y)$.
inline void bmm_geom2d_from_polar(double* const cart,
    double const* const polar) {
  cart[0] = polar[0] * cos(polar[1]);
  cart[1] = polar[0] * sin(polar[1]);
}

// The call `bmm_geom2d_from_polar2(cart, polar2)` maps
// the two-dimensional polar coordinates `polar2` of the form $(r^2, \phi)$
// to the corresponding cartesian coordinates `cart` of the form $(x, y)$.
inline void bmm_geom2d_from_polar2(double* const cart,
    double const* const polar2) {
  double const r = sqrt(polar2[0]);

  cart[0] = r * cos(polar2[1]);
  cart[1] = r * sin(polar2[1]);
}

// The call `bmm_geom2d_pdiff(rdiff, r0, r1, rper)` sets
// the two-dimensional vector `rdiff` to the `rper`-periodic difference
// between the vectors `r0` and `r1` by following the minimum image convention.
inline void bmm_geom2d_pdiff(double* const rdiff,
    double const* const r0, double const* const r1,
    double const* const rper) {
  for (size_t idim = 0; idim < 2; ++idim)
    rdiff[idim] = bmm_fp_swrap(r1[idim] - r0[idim], rper[idim]);
}

// The call `bmm_geom2d_pdist2(r0, r1, rper)` returns
// the `rper`-periodic distance $r^2$
// between the vectors `r0` and `r1` by following the minimum image convention.
__attribute__ ((__pure__))
inline double bmm_geom2d_pdist2(double const* const r0, double const* const r1,
    double const* const rper) {
  double r[2];
  bmm_geom2d_pdiff(r, r0, r1, rper);

  return bmm_geom2d_norm2(r);
}

// The call `bmm_geom2d_pdist(r0, r1, rper)` returns
// the `rper`-periodic distance $r$
// between the vectors `r0` and `r1` by following the minimum image convention.
__attribute__ ((__pure__))
inline double bmm_geom2d_pdist(double const* const r0, double const* const r1,
    double const* const rper) {
  double r[2];
  bmm_geom2d_pdiff(r, r0, r1, rper);

  return bmm_geom2d_norm(r);
}

// The call `bmm_geom2d_pangle(r0, r1, rper)` returns
// the `rper`-periodic signed angle $\phi$
// between the two-dimensional vectors `r0` and `r1`
// according to the minimum image convention.
__attribute__ ((__pure__))
inline double bmm_geom2d_pangle(double const* const r0, double const* const r1,
    double const* const rper) {
  double r[2];
  bmm_geom2d_pdiff(r, r0, r1, rper);

  return bmm_geom2d_angle(r);
}

#endif
