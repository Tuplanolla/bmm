/// Floating-point operations.

#ifndef BMM_FP_H
#define BMM_FP_H

#include <math.h>
#include <stddef.h>

#include "common.h"
#include "ext.h"

#ifndef M_PI_4
#define M_PI_4 0.7853981633974483
#endif

#ifndef M_PI_2
#define M_PI_2 1.5707963267948966
#endif

#ifndef M_PI
#define M_PI 3.141592653589793
#endif

#ifndef M_2PI
#define M_2PI 6.283185307179586
#endif

#ifndef M_4PI
#define M_4PI 12.566370614359172
#endif

/// The call `bmm_fp_lexcmp(x, y, n)`
/// performs `bmm_fp_cmp` on the arrays `x` and `y` of length `n`
/// in lexicographic order.
__attribute__ ((__pure__))
inline int bmm_fp_lexcmp(double const *restrict const x,
    double const *restrict const y, size_t const n) {
  for (size_t i = 0; i < n; ++i)
    switch ($(bmm_cmp, double)(x[i], y[i])) {
      case -1:
        return -1;
      case 1:
        return 1;
    }

  return 0;
}

/// The call `bmm_fp_rt(x, n)` returns the `n`th root of `x`.
/// This is analogous to `bmm_size_firt` or `bmm_size_cirt`.
__attribute__ ((__const__, __pure__))
inline double bmm_fp_rt(double const x, size_t const n) {
  return pow(x, 1.0 / (double) n);
}

/// The call `bmm_fp_log(x, y)` returns the the base `y` logarithm of `x`.
/// This is analogous to `bmm_size_flog` or `bmm_size_clog`.
__attribute__ ((__const__, __pure__))
inline double bmm_fp_log(double const x, double const y) {
  return log2(x) / log2(y);
}

/// The call `bmm_fp_clamp(x, a, b)` returns
///
/// * `x` if `a <= x < b`,
/// * `a` if `x < a` and
/// * `b` if `x >= b`.
__attribute__ ((__const__, __pure__))
inline double bmm_fp_clamp(double const x, double const a, double const b) {
  return x < a ? a : x >= b ? b : x;
}

/// The call `bmm_fp_sclamp(x, b)` returns
///
/// * `x` if `-b / 2 <= x < b / 2`,
/// * `-b / 2` if `x < b / 2` and
/// * `b` if `x >= b / 2`.
__attribute__ ((__const__, __pure__))
inline double bmm_fp_sclamp(double const x, double const b) {
  double const a = b / 2.0;

  return x < -a ? a : x >= a ? a : x;
}

/// The call `bmm_fp_uclamp(x, b)` returns
///
/// * `x` if `0 <= x < b`,
/// * `0` if `x < 0` and
/// * `b` if `x > b`.
///
/// This is analogous to `bmm_size_uclamp`.
__attribute__ ((__const__, __pure__))
inline double bmm_fp_uclamp(double const x, double const b) {
  return x < 0.0 ? 0.0: x >= b ? b : x;
}

/// The call `bmm_fp_min(x, n)`
/// returns the minimum of the array `x` of length `n`.
__attribute__ ((__pure__))
inline double bmm_fp_min(double const *const x, size_t const n) {
  double y = (double) INFINITY;

  for (size_t i = 0; i < n; ++i)
    y = fmin(y, x[i]);

  return y;
}

/// The call `bmm_fp_max(x, n)`
/// returns the maximum of the array `x` of length `n`.
__attribute__ ((__pure__))
inline double bmm_fp_max(double const *const x, size_t const n) {
  double y = -(double) INFINITY;

  for (size_t i = 0; i < n; ++i)
    y = fmax(y, x[i]);

  return y;
}

/// The call `y = bmm_fp_lerp(x, x0, x1, y0, y1)`
/// solves the linear interpolation equation
/// `(x - x0) / (x1 - x0) == (y - y0) / (y1 - y0)` for `y`.
__attribute__ ((__const__, __pure__))
inline double bmm_fp_lerp(double const x,
    double const x0, double const x1, double const y0, double const y1) {
  return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
}

/// The call `y = bmm_fp_lorp(x, x0, x1, y0, y1)`
/// solves the logarithmic interpolation equation and is equivalent to
/// `y = log(bmm_fp_lerp(exp(x), exp(x0), exp(x1), exp(y0), exp(y1)))`.
__attribute__ ((__const__, __pure__))
inline double bmm_fp_lorp(double const x,
    double const x0, double const x1, double const y0, double const y1) {
  return log(bmm_fp_lerp(exp(x), exp(x0), exp(x1), exp(y0), exp(y1)));
}

/// The call `bmm_fp_percent(x, y)`
/// returns the approximate percentage of `x` in `y`.
__attribute__ ((__const__, __pure__))
inline double bmm_fp_percent(double const x, double const y) {
  return y != 0.0 ? (x / y) * 100.0 : 100.0;
}

/// The call `n = bmm_fp_iclerp(x, x0, x1, n0, n1)`
/// is a discrete clamped version of `bmm_fp_lerp`.
__attribute__ ((__const__, __pure__))
inline size_t bmm_fp_iclerp(double const x,
    double const x0, double const x1, size_t const n0, size_t const n1) {
  dynamic_assert(n0 > 0, "Lower bound could wrap");
  // dynamic_assert(n1 <= SIZE_MAX, "Upper bound could wrap");

  double const y0 = (double) n0;
  double const y1 = (double) n1;
  double const y = bmm_fp_lerp(x, x0, x1, y0, y1);

  if (y < y0)
    return n0 - 1;
  else if (y >= y1)
    return n1;

  size_t const n = (size_t) y;
  dynamic_assert(n >= n0 && n < n1, "Invalid truncation");

  return n;
}

/// The call `bmm_fp_sumdist(x, a, b, c, d)`
/// calculates the value of the sum distribution
/// of two uniform distributions from `a` to `b` and from `c` to `d` at `x`.
__attribute__ ((__const__, __pure__))
inline double bmm_fp_sumdist(double const x,
    double const a, double const b,
    double const c, double const d) {
  double const ad = a + d;
  double const bc = b + c;

  if (ad < bc)
    return bmm_fp_sumdist(x, c, d, a, b);

  double const ac = a + c;
  double const bd = b + d;

  double const ba = b - a;
  double const dc = d - c;
  double const badc = ba * dc;

  return ac < x && x < bc ? (x - ac) / badc :
    bc < x && x < ad ? ba / badc :
    ad < x && x < bd ? -(x - bd) / badc :
    0.0;
}

/// The call `bmm_fp_proddist(x, a, b, c, d)`
/// calculates the value of the product distribution
/// of two uniform distributions from `a` to `b` and from `c` to `d` at `x`.
__attribute__ ((__const__, __pure__))
inline double bmm_fp_proddist(double const x,
    double const a, double const b,
    double const c, double const d) {
  double const ad = a * d;
  double const bc = b * c;

  if (ad < bc)
    return bmm_fp_proddist(x, c, d, a, b);

  double const ac = a * c;
  double const bd = b * d;

  double const ba = b - a;
  double const dc = d - c;
  double const badc = ba * dc;

  return ac < x && x < bc ? (log(x) - log(ac)) / badc :
    bc < x && x < ad ? (log(b) - log(a)) / badc :
    ad < x && x < bd ? log(bd / x) / badc :
    0.0;
}


#endif
