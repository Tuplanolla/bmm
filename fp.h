/// Floating-point operations.
#ifndef BMM_FP_H
#define BMM_FP_H

#include <math.h>

#include "ext.h"

#ifndef M_PI
#define M_PI 3.141592653589793
#endif

#ifndef M_2PI
#define M_2PI 6.283185307179586
#endif

/// This structure holds the quotient and remainder of a division
/// in unspecified order.
typedef struct {
  double quot;
  double rem;
} bmm_fp_div_t;

/// The call `z = bmm_fp_div(x, y)`
/// solves the division equation `z.quot * y + z.rem == x` for `z`,
/// where `z.quot` is the quotient and `z.rem` is the remainder
/// of the expression `x / y`.
/// This is analogous to `div` or `bmm_size_div`.
__attribute__ ((__const__, __pure__))
inline bmm_fp_div_t bmm_fp_div(double const x, double const y) {
  bmm_fp_div_t const qr = {
    .quot = trunc(x / y),
    .rem = fmod(x, y)
  };

  return qr;
}

/// The call `bmm_fp_cmp(x, y)` returns
///
/// * `-1` if `x < y`,
/// * `1` if `x > y` and
/// * `0` otherwise.
///
/// This is analogous to `bmm_size_cmp`.
__attribute__ ((__const__, __pure__))
inline int bmm_fp_cmp(double const x, double const y) {
  return x < y ? -1 : x > y ? 1 : 0;
}

/// The call `bmm_fp_identity(x)` returns `x`.
/// This is analogous to `bmm_size_identity`.
__attribute__ ((__const__, __pure__))
inline double bmm_fp_identity(double const x) {
  return x;
}

/// The call `bmm_fp_constant(x, y)` returns `x`.
/// This is analogous to `bmm_size_constant`.
__attribute__ ((__const__, __pure__))
inline double bmm_fp_constant(double const x,
    __attribute__ ((__unused__)) double const y) {
  return x;
}

/// The call `bmm_fp_zero(x)` returns `0`.
/// This is analogous to `bmm_size_zero`.
__attribute__ ((__const__, __pure__))
inline double bmm_fp_zero(__attribute__ ((__unused__)) double const x) {
  return 0.0;
}

/// The call `bmm_fp_one(x)` returns `1`.
/// This is analogous to `bmm_size_one`.
__attribute__ ((__const__, __pure__))
inline double bmm_fp_one(__attribute__ ((__unused__)) double const x) {
  return 1.0;
}

/// The call `bmm_fp_midpoint(x, y)`
/// returns the arithmetic mean of `x` and `y`.
/// This is analogous to `bmm_size_midpoint`.
__attribute__ ((__const__, __pure__))
inline double bmm_fp_midpoint(double const x, double const y) {
  return (x + y) / 2.0;
}

/// The call `bmm_fp_sq(x)` returns `x` squared.
__attribute__ ((__const__, __pure__))
inline double bmm_fp_sq(double const x) {
  return x * x;
}

/// The call `bmm_fp_rt(x, y)` returns the `y`th root of `x`.
/// This is analogous to `bmm_size_firt` or `bmm_size_cirt`.
__attribute__ ((__const__, __pure__))
inline double bmm_fp_rt(double const x, double const y) {
  return pow(x, 1.0 / y);
}

/// The call `bmm_fp_log(x, y)` returns the the base `y` logarithm of `x`.
/// This is analogous to `bmm_size_filog` or `bmm_size_cilog`.
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

/// The call `z = bmm_fp_wrap(x, a, b)`
/// solves the periodic equation `z == x - a + n * a` for `z`,
/// where `a <= z < b` and `n` is some integer.
__attribute__ ((__const__, __pure__))
inline double bmm_fp_wrap(double const x, double const a, double const b) {
  double const c = b - a;

  return x - c * floor((x - a) / c);
}

/// The call `z = bmm_fp_swrap(x, b)`
/// solves the periodic equation `z == x + n * b` for `z`,
/// where `-b / 2 <= z < b / 2` and `n` is some integer.
/// The `s` prefix means signed or symmetric.
__attribute__ ((__const__, __pure__))
inline double bmm_fp_swrap(double const x, double const b) {
  return x - b * nearbyint(x / b);
}

/// The call `z = bmm_fp_uwrap(x, b)`
/// solves the periodic equation `z == x + n * b` for `z`,
/// where `0 <= z < b` and `n` is some integer.
/// This is analogous to `bmm_size_uwrap`.
/// The `u` prefix means unsigned or unsymmetric (asymmetric).
__attribute__ ((__const__, __pure__))
inline double bmm_fp_uwrap(double const x, double const b) {
  return x - b * floor(x / b);
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

#endif
