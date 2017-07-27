/// Floating-point operations.

#ifndef BMM_FP_H
#define BMM_FP_H

#include <math.h>
#include <stddef.h>

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
/// This is analogous to `div` or `type(bmm_div, size_t)`.
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
/// This is analogous to `type(bmm_cmp, size_t)`.
__attribute__ ((__const__, __pure__))
inline int bmm_fp_cmp(double const x, double const y) {
  return x < y ? -1 : x > y ? 1 : 0;
}

/// The call `bmm_fp_lexcmp(x, y, n)`
/// performs `bmm_fp_cmp` on the arrays `x` and `y` of length `n`
/// in lexicographic order.
__attribute__ ((__pure__))
inline int bmm_fp_lexcmp(double const *restrict const x,
    double const *restrict const y, size_t const n) {
  for (size_t i = 0; i < n; ++i)
    switch (bmm_fp_cmp(x[i], y[i])) {
      case -1:
        return -1;
      case 1:
        return 1;
    }

  return 0;
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

/// The call `bmm_fp_cb(x)` returns `x` cubed.
__attribute__ ((__const__, __pure__))
inline double bmm_fp_cb(double const x) {
  return x * x * x;
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

/// The call `z = bmm_fp_wrap(x, a, b)`
/// solves the periodic equation `z == x - a + k * a` for `z`,
/// where `a <= z < b` and `k` is some integer.
__attribute__ ((__const__, __pure__))
inline double bmm_fp_wrap(double const x, double const a, double const b) {
  double const c = b - a;

  return x - c * floor((x - a) / c);
}

/// The call `z = bmm_fp_swrap(x, b)`
/// solves the periodic equation `z == x + k * b` for `z`,
/// where `-b / 2 <= z < b / 2` and `k` is some integer.
/// The `s` prefix means signed or symmetric.
__attribute__ ((__const__, __pure__))
inline double bmm_fp_swrap(double const x, double const b) {
  return x - b * nearbyint(x / b);
}

/// The call `z = bmm_fp_uwrap(x, b)`
/// solves the periodic equation `z == x + k * b` for `z`,
/// where `0 <= z < b` and `k` is some integer.
/// This is analogous to `type(bmm_uwrap, size_t)`.
/// The `u` prefix means unsigned or unsymmetric (asymmetric).
__attribute__ ((__const__, __pure__))
inline double bmm_fp_uwrap(double const x, double const b) {
  return x - b * floor(x / b);
}

/// The call `bmm_fp_sum(x, n)`
/// returns the sum of the array `x` of length `n`.
__attribute__ ((__pure__))
inline double bmm_fp_sum(double const *const x, size_t const n) {
  double y = 0.0;

  for (size_t i = 0; i < n; ++i)
    y += x[i];

  return y;
}

/// The call `bmm_fp_prod(x, n)`
/// returns the product of the array `x` of length `n`.
__attribute__ ((__pure__))
inline double bmm_fp_prod(double const *const x, size_t const n) {
  double y = 1.0;

  for (size_t i = 0; i < n; ++i)
    y *= x[i];

  return y;
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

/// The call `bmm_fp_lfold(f, x, n, z, ptr)`
/// folds the procedure `f` over the array `x` of length `n`.
/// by starting from the left with `z`.
__attribute__ ((__nonnull__ (1, 2)))
inline double bmm_fp_lfold(double (*const f)(double, double, void *),
    double const *restrict const x, size_t const n,
    double z, void *restrict const ptr) {
  for (size_t i = 0; i < n; ++i)
    z = f(x[i], z, ptr);

  return z;
}

/// The call `bmm_fp_rfold(f, x, n, z, ptr)`
/// folds the procedure `f` over the array `x` of length `n`.
/// by starting from the right with `z`.
__attribute__ ((__nonnull__ (1, 2)))
inline double bmm_fp_rfold(double (*const f)(double, double, void *),
    double const *restrict const x, size_t const n,
    double z, void *restrict const ptr) {
  for (size_t i = 0; i < n; ++i)
    z = f(x[n - 1 - i], z, ptr);

  return z;
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

/// The call `bmm_fp_pow(x, n)` returns the `n`th power of `x`.
/// This is analogous to `pow`.
__attribute__ ((__const__, __pure__))
inline double bmm_fp_pow(double const x, size_t const n) {
  return pow(x, (double) n);
}

/// The call `bmm_fp_fact(n, k)`
/// returns the `k`-factorial of `n`.
__attribute__ ((__const__, __pure__))
inline double bmm_fp_fact(size_t const n, size_t const k) {
  double x = 1.0;

  for (size_t i = n; i > 1; i -= k)
    x *= (double) i;

  return x;
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

#endif
