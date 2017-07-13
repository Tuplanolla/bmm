#ifndef BMM_GEOM2D_H
/// Planar geometry.
#define BMM_GEOM2D_H

#include <math.h>
#include <stdbool.h>
#include <stddef.h>

#include "ext.h"
#include "fp.h"
#include "ival.h"

// TODO Add __nonnull__.

/// The call `bmm_geom2d_dot(x)`
/// returns the dot product of the vectors `x` and `y`.
__attribute__ ((__pure__))
inline double bmm_geom2d_dot(double const* restrict const x,
    double const* restrict const y) {
  double d = 0.0;
  for (size_t idim = 0; idim < 2; ++idim)
    d += x[idim] * y[idim];

  return d;
}

/// The call `bmm_geom2d_norm2(x)`
/// returns the length $r^2$ of the vector `x`.
__attribute__ ((__pure__))
inline double bmm_geom2d_norm2(double const* const x) {
  double d = 0.0;
  for (size_t idim = 0; idim < 2; ++idim)
    d += bmm_fp_sq(x[idim]);

  return d;
}

/// The call `bmm_geom2d_norm(x)`
/// returns the length $r$ of the vector `x`.
__attribute__ ((__pure__))
inline double bmm_geom2d_norm(double const* const x) {
  return sqrt(bmm_geom2d_norm2(x));
}

/// The call `bmm_geom2d_dir(x)`
/// returns the signed angle $\\phi$ of the vector `x`.
__attribute__ ((__pure__))
inline double bmm_geom2d_dir(double const* const x) {
  return atan2(x[1], x[0]);
}

/// The call `bmm_geom2d_redir(phi)`
/// returns the reversed signed angle $\\phi + \\twopi / 2$ of the angle `phi`.
__attribute__ ((__const__, __pure__))
inline double bmm_geom2d_redir(double const phi) {
  return phi + M_2PI / 2;
}

/// The call `bmm_geom2d_add(y, x0, x1)`
/// sets the vector `y` to the vector `x0` plus `x1`.
inline void bmm_geom2d_add(double* restrict const y,
    double const* restrict const x0, double const* restrict const x1) {
  for (size_t idim = 0; idim < 2; ++idim)
    y[idim] = x0[idim] + x1[idim];
}

/// The call `bmm_geom2d_addto(y, x)`
/// sets the vector `y` to the vector `x` plus `y`.
inline void bmm_geom2d_addto(double* restrict const y,
    double const* restrict const x) {
  for (size_t idim = 0; idim < 2; ++idim)
    y[idim] += x[idim];
}

/// The call `bmm_geom2d_scaler(y, x, a)`
/// sets the vector `y` to the vector `x` scaled by `a`.
inline void bmm_geom2d_scale(double* restrict const y,
    double const* restrict const x, double const a) {
  for (size_t idim = 0; idim < 2; ++idim)
    y[idim] = a * x[idim];
}

/// The call `bmm_geom2d_normal(y, x)`
/// finds the normal vector `y` of the vector `x`.
inline void bmm_geom2d_normal(double* restrict const y,
    double const* restrict const x) {
  bmm_geom2d_scale(y, x, bmm_geom2d_norm(x));
}

/// The call `bmm_geom2d_project(z, x, y)`
/// finds the projection `z` of the vector `x` on the normal vector `y`.
inline void bmm_geom2d_project(double* restrict const z,
    double const* restrict const x, double const* restrict const y) {
  bmm_geom2d_scale(z, y, bmm_geom2d_dot(x, y));
}

/// The call `bmm_geom2d_rperp(y, x)`
/// finds the right-handed perpendicular vector `y` of the form $(-y, x)$
/// given the vector `x` of the form $(x, y)$.
inline void bmm_geom2d_rperp(double* restrict const y,
    double const* restrict const x) {
  y[0] = -x[1];
  y[1] = x[0];
}

/// The call `bmm_geom2d_lperp(y, x)`
/// finds the left-handed perpendicular vector `y` of the form $(y, -x)$
/// given the vector `x` of the form $(x, y)$.
inline void bmm_geom2d_lperp(double* restrict const y,
    double const* restrict const x) {
  y[1] = -x[0];
  y[0] = x[1];
}

/// The call `bmm_geom2d_to_polar2(polar2, cart)`
/// maps the cartesian coordinates `cart` of the form $(x, y)$
/// to the corresponding polar coordinates `polar2` of the form $(r^2, \\phi)$.
inline void bmm_geom2d_to_polar2(double* restrict const polar2,
    double const* restrict const cart) {
  polar2[0] = bmm_geom2d_norm2(cart);
  polar2[1] = bmm_geom2d_dir(cart);
}

/// The call `bmm_geom2d_to_polar(polar, cart)`
/// maps the cartesian coordinates `cart` of the form $(x, y)$
/// to the corresponding polar coordinates `polar` of the form $(r, \\phi)$.
inline void bmm_geom2d_to_polar(double* restrict const polar,
    double const* restrict const cart) {
  polar[0] = bmm_geom2d_norm(cart);
  polar[1] = bmm_geom2d_dir(cart);
}

/// The call `bmm_geom2d_from_polar(cart, polar)`
/// maps the polar coordinates `polar` of the form $(r, \\phi)$
/// to the corresponding cartesian coordinates `cart` of the form $(x, y)$.
inline void bmm_geom2d_from_polar(double* restrict const cart,
    double const* restrict const polar) {
  cart[0] = polar[0] * cos(polar[1]);
  cart[1] = polar[0] * sin(polar[1]);
}

/// The call `bmm_geom2d_from_polar2(cart, polar2)`
/// maps the polar coordinates `polar2` of the form $(r^2, \\phi)$
/// to the corresponding cartesian coordinates `cart` of the form $(x, y)$.
inline void bmm_geom2d_from_polar2(double* restrict const cart,
    double const* restrict const polar2) {
  double const r = sqrt(polar2[0]);

  cart[0] = r * cos(polar2[1]);
  cart[1] = r * sin(polar2[1]);
}

/// The call `bmm_geom2d_diff(xdiff, x0, x1)`
/// sets the vector `xdiff` to the difference
/// between the vectors `x0` and `x1`.
inline void bmm_geom2d_diff(double* restrict const xdiff,
    double const* restrict const x0, double const* restrict const x1) {
  for (size_t idim = 0; idim < 2; ++idim)
    xdiff[idim] = x1[idim] - x0[idim];
}

/// The call `bmm_geom2d_dist2(x0, x1)`
/// returns the distance $r^2$ between the vectors `x0` and `x1`.
__attribute__ ((__pure__))
inline double bmm_geom2d_dist2(double const* restrict const x0,
    double const* restrict const x1) {
  double x[2];
  bmm_geom2d_diff(x, x0, x1);

  return bmm_geom2d_norm2(x);
}

/// The call `bmm_geom2d_dist(x0, x1)`
/// returns the distance $r$ between the vectors `x0` and `x1`.
__attribute__ ((__pure__))
inline double bmm_geom2d_dist(double const* restrict const x0,
    double const* restrict const x1) {
  double x[2];
  bmm_geom2d_diff(x, x0, x1);

  return bmm_geom2d_norm(x);
}

/// The call `bmm_geom2d_angle(x0, x1)`
/// returns the signed angle $\\phi$
/// between the vectors `x0` and `x1`
/// according to the minimum image convention.
__attribute__ ((__pure__))
inline double bmm_geom2d_angle(double const* restrict const x0,
    double const* restrict const x1) {
  double x[2];
  bmm_geom2d_diff(x, x0, x1);

  return bmm_geom2d_dir(x);
}

/// The call `bmm_geom2d_pdiff(xdiff, x0, x1, xper)`
/// sets the vector `xdiff` to the `xper`-periodic difference
/// between the vectors `x0` and `x1`
/// by following the minimum image convention.
inline void bmm_geom2d_pdiff(double* restrict const xdiff,
    double const* restrict const x0, double const* restrict const x1,
    double const* restrict const xper) {
  bmm_geom2d_diff(xdiff, x0, x1);

  for (size_t idim = 0; idim < 2; ++idim)
    xdiff[idim] = bmm_fp_swrap(xdiff[idim], xper[idim]);
}

/// The call `bmm_geom2d_pdist2(x0, x1, xper)`
/// returns the `xper`-periodic distance $r^2$
/// between the vectors `x0` and `x1`
/// by following the minimum image convention.
__attribute__ ((__pure__))
inline double bmm_geom2d_pdist2(double const* restrict const x0,
    double const* restrict const x1, double const* restrict const xper) {
  double x[2];
  bmm_geom2d_pdiff(x, x0, x1, xper);

  return bmm_geom2d_norm2(x);
}

/// The call `bmm_geom2d_pdist(x0, x1, xper)`
/// returns the `xper`-periodic distance $r$
/// between the vectors `x0` and `x1`
/// by following the minimum image convention.
__attribute__ ((__pure__))
inline double bmm_geom2d_pdist(double const* restrict const x0,
    double const* restrict const x1, double const* restrict const xper) {
  double x[2];
  bmm_geom2d_pdiff(x, x0, x1, xper);

  return bmm_geom2d_norm(x);
}

/// The call `bmm_geom2d_pangle(x0, x1, xper)`
/// returns the `xper`-periodic signed angle $\\phi$
/// between the vectors `x0` and `x1`
/// according to the minimum image convention.
__attribute__ ((__pure__))
inline double bmm_geom2d_pangle(double const* restrict const x0,
    double const* restrict const x1, double const* restrict const xper) {
  double x[2];
  bmm_geom2d_pdiff(x, x0, x1, xper);

  return bmm_geom2d_dir(x);
}

/// The call `bmm_geom2d_cpdiff(xdiff, x0, x1, xper, per)`
/// sets the vector `xdiff` to the `per`-conditional `xper`-periodic difference
/// between the vectors `x0` and `x1`
/// by following the minimum image convention.
inline void bmm_geom2d_cpdiff(double* restrict const xdiff,
    double const* restrict const x0, double const* restrict const x1,
    double const* restrict const xper, bool const* restrict const per) {
  bmm_geom2d_diff(xdiff, x0, x1);

  for (size_t idim = 0; idim < 2; ++idim)
    if (per[idim])
      xdiff[idim] = bmm_fp_swrap(xdiff[idim], xper[idim]);
}

/// The call `bmm_geom2d_cpdist2(x0, x1, xper, per)`
/// returns the `per`-conditional `xper`-periodic distance $r^2$
/// between the vectors `x0` and `x1`
/// by following the minimum image convention.
__attribute__ ((__pure__))
inline double bmm_geom2d_cpdist2(double const* restrict const x0,
    double const* restrict const x1, double const* restrict const xper,
    bool const* restrict const per) {
  double x[2];
  bmm_geom2d_cpdiff(x, x0, x1, xper, per);

  return bmm_geom2d_norm2(x);
}

/// The call `bmm_geom2d_cpdist(x0, x1, xper, per)`
/// returns the `per`-conditional `xper`-periodic distance $r$
/// between the vectors `x0` and `x1`
/// by following the minimum image convention.
__attribute__ ((__pure__))
inline double bmm_geom2d_cpdist(double const* restrict const x0,
    double const* restrict const x1, double const* restrict const xper,
    bool const* restrict const per) {
  double x[2];
  bmm_geom2d_cpdiff(x, x0, x1, xper, per);

  return bmm_geom2d_norm(x);
}

/// The call `bmm_geom2d_cpangle(x0, x1, xper, per)`
/// returns the `per`-conditional `xper`-periodic signed angle $\\phi$
/// between the vectors `x0` and `x1`
/// according to the minimum image convention.
__attribute__ ((__pure__))
inline double bmm_geom2d_cpangle(double const* restrict const x0,
    double const* restrict const x1, double const* restrict const xper,
    bool const* restrict const per) {
  double x[2];
  bmm_geom2d_cpdiff(x, x0, x1, xper, per);

  return bmm_geom2d_dir(x);
}

__attribute__ ((__nonnull__, __pure__))
inline double bmm_geom2d_quadal(double const* restrict const x,
    double const r,
    double const* restrict const xper, bool const* const per) {
  if (r <= 0.0)
    return 0.0;

  double y[2];
  bmm_geom2d_diff(y, x, xper);

  if (bmm_geom2d_norm2(y) <= bmm_fp_pow(r, 2))
    return 0.0;

  double a[2];
  a[0] = per[0] || y[0] >= r ? 0.0 : acos(y[0] / r);
  a[1] = per[1] || y[1] >= r ? M_PI_2 : asin(y[1] / r);

  return r * bmm_ival_length(a);
}

__attribute__ ((__nonnull__, __pure__))
inline double bmm_geom2d_shellal(double const* restrict const x,
    double const r,
    double const* restrict const xper, bool const* const per) {
  // TODO Inside?
  double al = 0.0;

  double y[2];

  // TODO Reflect explicitly?
  y[0] = x[0];
  y[1] = x[1];
  al += bmm_geom2d_quadal(y, r, xper, per);

  y[0] = xper[0] - x[0];
  y[1] = x[1];
  al += bmm_geom2d_quadal(y, r, xper, per);

  y[0] = x[0];
  y[1] = xper[1] - x[1];
  al += bmm_geom2d_quadal(y, r, xper, per);

  y[0] = xper[0] - x[0];
  y[1] = xper[1] - x[1];
  al += bmm_geom2d_quadal(y, r, xper, per);

  return al;
}

#endif
