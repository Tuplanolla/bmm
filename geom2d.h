// Planar geometry.
#ifndef BMM_GEOM2D_H
#define BMM_GEOM2D_H

#include "ext.h"
#include "fp.h"
#include <math.h>
#include <stddef.h>

// The call `bmm_geom2d_norm2(r)` returns the length $r^2$
// of the vector `r`.
__attribute__ ((__pure__))
inline double bmm_geom2d_norm2(double const* const r) {
  double d = 0.0;
  for (size_t idim = 0; idim < 2; ++idim)
    d += bmm_fp_sq(r[idim]);

  return d;
}

// The call `bmm_geom2d_norm(r)` returns the length $r$
// of the vector `r`.
__attribute__ ((__pure__))
inline double bmm_geom2d_norm(double const* const r) {
  return sqrt(bmm_geom2d_norm2(r));
}

// The call `bmm_geom2d_dir(r)` returns the signed angle $\phi$
// of the vector `r`.
__attribute__ ((__pure__))
inline double bmm_geom2d_dir(double const* const r) {
  return atan2(r[1], r[0]);
}

// The call `bmm_geom2d_add(a, r0, r1)`
// sets the vector `a` to the vector `r0` plus `r1`.
inline void bmm_geom2d_add(double* const a,
    double const* const r0, double const* const r1) {
  for (size_t idim = 0; idim < 2; ++idim)
    a[idim] = r0[idim] + r1[idim];
}

// The call `bmm_geom2d_scale(s, r, x)`
// sets the vector `s` to the vector `r` scaled by `x`.
inline void bmm_geom2d_scale(double* const s,
    double const* const r, double const x) {
  for (size_t idim = 0; idim < 2; ++idim)
    s[idim] = x * r[idim];
}

// The call `bmm_geom2d_normal(n, r)` finds
// the normal vector `n` of the vector `r`.
inline void bmm_geom2d_normal(double* const n, double const* const r) {
  bmm_geom2d_scale(n, r, bmm_geom2d_norm(r));
}

// The call `bmm_geom2d_rperpr(p, r)` finds
// the right-handed perpendicular vector `p` of the form $(-y, x)$
// given the vector `r` of the form $(x, y)$.
inline void bmm_geom2d_rperpr(double* restrict const p,
    double const* restrict const r) {
  p[0] = -r[1];
  p[1] = r[0];
}

// The call `bmm_geom2d_rperp(p, r)` finds
// the right-handed perpendicular vector `p` of the form $(-y, x)$
// given the vector `r` of the form $(x, y)$.
inline void bmm_geom2d_rperp(double* const p, double const* const r) {
  double const x = r[0];
  p[0] = -r[1];
  p[1] = x;
}

// The call `bmm_geom2d_lperpr(p, r)` finds
// the left-handed perpendicular vector `p` of the form $(y, -x)$
// given the vector `r` of the form $(x, y)$.
inline void bmm_geom2d_lperpr(double* restrict const p,
    double const* restrict const r) {
  p[1] = -r[0];
  p[0] = r[1];
}

// The call `bmm_geom2d_lperp(p, r)` finds
// the left-handed perpendicular vector `p` of the form $(y, -x)$
// given the vector `r` of the form $(x, y)$.
inline void bmm_geom2d_lperp(double* const p, double const* const r) {
  double const y = r[1];
  p[1] = -r[0];
  p[0] = y;
}

// The call `bmm_geom2d_to_polar2(polar2, cart)` maps
// the cartesian coordinates `cart` of the form $(x, y)$
// to the corresponding polar coordinates `polar2` of the form $(r^2, \phi)$.
inline void bmm_geom2d_to_polar2(double* const polar2,
    double const* const cart) {
  polar2[0] = bmm_geom2d_norm2(cart);
  polar2[1] = bmm_geom2d_dir(cart);
}

// The call `bmm_geom2d_to_polar(polar, cart)` maps
// the cartesian coordinates `cart` of the form $(x, y)$
// to the corresponding polar coordinates `polar` of the form $(r, \phi)$.
inline void bmm_geom2d_to_polar(double* const polar,
    double const* const cart) {
  polar[0] = bmm_geom2d_norm(cart);
  polar[1] = bmm_geom2d_dir(cart);
}

// The call `bmm_geom2d_from_polar(cart, polar)` maps
// the polar coordinates `polar` of the form $(r, \phi)$
// to the corresponding cartesian coordinates `cart` of the form $(x, y)$.
inline void bmm_geom2d_from_polar(double* const cart,
    double const* const polar) {
  cart[0] = polar[0] * cos(polar[1]);
  cart[1] = polar[0] * sin(polar[1]);
}

// The call `bmm_geom2d_from_polar2(cart, polar2)` maps
// the polar coordinates `polar2` of the form $(r^2, \phi)$
// to the corresponding cartesian coordinates `cart` of the form $(x, y)$.
inline void bmm_geom2d_from_polar2(double* const cart,
    double const* const polar2) {
  double const r = sqrt(polar2[0]);

  cart[0] = r * cos(polar2[1]);
  cart[1] = r * sin(polar2[1]);
}

// The call `bmm_geom2d_diff(rdiff, r0, r1)`
// sets the vector `rdiff` to the difference between the vectors `r0` and `r1`.
inline void bmm_geom2d_diff(double* const rdiff,
    double const* const r0, double const* const r1) {
  for (size_t idim = 0; idim < 2; ++idim)
    rdiff[idim] = r1[idim] - r0[idim];
}

// The call `bmm_geom2d_dist2(r0, r1)`
// returns the distance $r^2$ between the vectors `r0` and `r1`.
__attribute__ ((__pure__))
inline double bmm_geom2d_dist2(double const* const r0,
    double const* const r1) {
  double r[2];
  bmm_geom2d_diff(r, r0, r1);

  return bmm_geom2d_norm2(r);
}

// The call `bmm_geom2d_dist(r0, r1)`
// returns the distance $r$ between the vectors `r0` and `r1`.
__attribute__ ((__pure__))
inline double bmm_geom2d_dist(double const* const r0, double const* const r1) {
  double r[2];
  bmm_geom2d_diff(r, r0, r1);

  return bmm_geom2d_norm(r);
}

// The call `bmm_geom2d_angle(r0, r1)`
// returns the signed angle $\phi$
// between the vectors `r0` and `r1`
// according to the minimum image convention.
__attribute__ ((__pure__))
inline double bmm_geom2d_angle(double const* const r0,
    double const* const r1) {
  double r[2];
  bmm_geom2d_diff(r, r0, r1);

  return bmm_geom2d_dir(r);
}

// The call `bmm_geom2d_pdiff(rdiff, r0, r1, rper)`
// sets the vector `rdiff` to the `rper`-periodic difference
// between the vectors `r0` and `r1` by following the minimum image convention.
inline void bmm_geom2d_pdiff(double* const rdiff,
    double const* const r0, double const* const r1,
    double const* const rper) {
  for (size_t idim = 0; idim < 2; ++idim)
    rdiff[idim] = bmm_fp_swrap(r1[idim] - r0[idim], rper[idim]);
}

// The call `bmm_geom2d_pdist2(r0, r1, rper)`
// returns the `rper`-periodic distance $r^2$
// between the vectors `r0` and `r1` by following the minimum image convention.
__attribute__ ((__pure__))
inline double bmm_geom2d_pdist2(double const* const r0, double const* const r1,
    double const* const rper) {
  double r[2];
  bmm_geom2d_pdiff(r, r0, r1, rper);

  return bmm_geom2d_norm2(r);
}

// The call `bmm_geom2d_pdist(r0, r1, rper)`
// returns the `rper`-periodic distance $r$
// between the vectors `r0` and `r1` by following the minimum image convention.
__attribute__ ((__pure__))
inline double bmm_geom2d_pdist(double const* const r0, double const* const r1,
    double const* const rper) {
  double r[2];
  bmm_geom2d_pdiff(r, r0, r1, rper);

  return bmm_geom2d_norm(r);
}

// The call `bmm_geom2d_pangle(r0, r1, rper)`
// returns the `rper`-periodic signed angle $\phi$
// between the vectors `r0` and `r1`
// according to the minimum image convention.
__attribute__ ((__pure__))
inline double bmm_geom2d_pangle(double const* const r0, double const* const r1,
    double const* const rper) {
  double r[2];
  bmm_geom2d_pdiff(r, r0, r1, rper);

  return bmm_geom2d_dir(r);
}

#endif
