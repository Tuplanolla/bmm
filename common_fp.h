/// Common operations for floating-point types.

#include "ext.h"
#include "wrap.h"

/// The call `bmm_quotrem(x, y)`
/// returns the quotient and remainder of `x` divided by `y`
/// in `qr` such that `qr.quot * y + qr.rem == x` and `qr.rem >= 0`.
/// If `y == 0` or `x` and `y` are infinite or `x` or `y` are not numbers,
/// the behavior is undefined.
__attribute__ ((__const__, __pure__))
inline $(bmm_quotrem_t, A) $(bmm_quotrem, A)(A const x, A const y) {
  A const q = $(trunc, A)(x / y);
  A const r = $(fmod, A)(x, y);
  A const s = r >= 0 ? 0 : y < 0 ? -1 : 1;

  return ($(bmm_quotrem_t, A)) {.quot = q - s, .rem = r + s * y};
}

/// The call `bmm_abs(x)`
/// returns the absolute value of `x`.
/// If `x` is not a number, the behavior is undefined.
__attribute__ ((__const__, __pure__))
inline A $(bmm_abs, A)(A const x) {
  return x < 0 ? -x : x;
}

/// The call `bmm_wrap(x, a, b)`
/// finds such `y` that `a <= y < b`
/// by shifting `x` by the appropriate number of multiples of `b - a`.
/// If `b <= a` or `x` is infinite or `x`, `a` or `b` are not numbers,
/// the behavior is undefined.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A $(bmm_wrap, A)(A const x, A const a, A const b) {
  dynamic_assert(b > a, "Invalid argument");

  A const c = b - a;

  return x - c * floor((x - a) / c);

  // The following implementation is more consistent with integers,
  // but potentially slower.
  // return $(bmm_rem, A)(x - a, b - a) + a;
}

/// The call `bmm_uwrap(x, b)`
/// finds such `y` that `0 <= y < b`
/// by shifting `x` by the appropriate number of multiples of `b`.
/// If `b <= 0` or `x` is infinite or `x` or `b` are not numbers,
/// the behavior is undefined.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A $(bmm_uwrap, A)(A const x, A const b) {
  dynamic_assert(b > 0, "Invalid argument");

  return x - b * floor(x / b);

  // The following implementation is more consistent with integers,
  // but potentially slower.
  // return $(bmm_rem, A)(x, b);
}

/// The call `bmm_swrap(x, c)`
/// finds such `y` that `-c / 2 <= y < c / 2`
/// by shifting `x` by the appropriate number of multiples of `c`.
/// If `c <= 0` or `x` is infinite or `x` or `c` are not numbers,
/// the behavior is undefined.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A $(bmm_swrap, A)(A const x, A const c) {
  dynamic_assert(c > 0, "Invalid argument");

  return x - c * floor(x / c + 0.5);

  // The following implementation is more consistent with integers,
  // but potentially slower.
  // A const d = c / 2;
  //
  // return $(bmm_wrap, A)(x, -d, d);
}

/// The call `bmm_uclamp(x, b)`
/// finds such `y` that `0 <= y <= b`
/// by shifting `x` by the smallest possible amount.
/// If `b < 0` or `x` is infinite or `x` or `b` are not numbers,
/// the behavior is undefined.
/// Overflows are impossible both internally and externally.
__attribute__ ((__const__, __pure__))
inline A $(bmm_uclamp, A)(A const x, A const b) {
  dynamic_assert(b >= 0, "Invalid argument");

  return x < 0 ? 0 : x > b ? b : x;
}

/// The call `bmm_sclamp(x, c)`
/// finds such `y` that `-c / 2 <= y <= c / 2`
/// by shifting `x` by the smallest possible amount.
/// If `c < 0` or `x` is infinite or `x` or `c` are not numbers,
/// the behavior is undefined.
__attribute__ ((__const__, __pure__))
inline A $(bmm_sclamp, A)(A const x, A const c) {
  dynamic_assert(c >= 0, "Invalid argument");

  A const b = c / 2;
  A const a = -b;

  return x < a ? a : x > b ? b : x;
}

/// The call `bmm_fact(x)`
/// returns the factorial of `x`.
/// If `x < 0` or `x` is infinite or `x` is not a number,
/// the behavior is undefined.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A $(bmm_fact, A)(A const x) {
  dynamic_assert(x >= 0, "Invalid argument");

  return tgamma(x + 1);
}

/// The call `bmm_resum2(x, y)`
/// returns the reciprocal sum of `x` and `y`.
/// If `x <= 0` or `y <= 0`, `x` or `y` are infinite or
/// `x` or `y` are not numbers, the behavior is undefined.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A $(bmm_resum2, A)(A const x, A const y) {
  dynamic_assert(x > 0, "Invalid argument");
  dynamic_assert(y > 0, "Invalid argument");

  return (x * y) / (x + y);

  // The following implementation is closer to the original definition,
  // but slower and less stable.
  // return 1 / (1 / x + 1 / y);
}

/// The call `bmm_pmean2(x, y, e)`
/// returns the power mean of `x` and `y` with the exponent `e`.
/// If `x` or `y` are infinite or
/// `x` or `y` are not numbers, the behavior is undefined.
__attribute__ ((__const__, __pure__))
inline A $(bmm_pmean2, A)(A const x, A const y, A const e) {
  A const z = ($(pow, A)(x, e) + $(pow, A)(y, e)) / 2;

  return $(pow, A)(z, 1 / e);
}

/// The call `bmm_amean2(x, y)`
/// returns the arithmetic mean of `x` and `y`.
/// It is equivalent to `bmm_pmean2(x, y, 1)`.
/// If `x` or `y` are infinite or
/// `x` or `y` are not numbers, the behavior is undefined.
__attribute__ ((__const__, __pure__))
inline A $(bmm_amean2, A)(A const x, A const y) {
  return (x + y) / 2;
}

/// The call `bmm_gmean2(x, y)`
/// returns the geometric mean of `x` and `y`.
/// It is equivalent to `bmm_pmean2(x, y, 0)` at the limit.
/// If `x < 0` or `y < 0`, `x` or `y` are infinite or
/// `x` or `y` are not numbers, the behavior is undefined.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A $(bmm_gmean2, A)(A const x, A const y) {
  dynamic_assert(x >= 0, "Invalid argument");
  dynamic_assert(y >= 0, "Invalid argument");

  return sqrt(x * y);
}

/// The call `bmm_hmean2(x, y)`
/// returns the harmonic mean of `x` and `y`.
/// It is equivalent to `bmm_pmean2(x, y, -1)`.
/// If `x <= 0` or `y <= 0`, `x` or `y` are infinite or
/// `x` or `y` are not numbers, the behavior is undefined.
__attribute__ ((__const__, __pure__))
inline A $(bmm_hmean2, A)(A const x, A const y) {
  return 2 * $(bmm_resum2, A)(x, y);
}
