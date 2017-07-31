/// Common operations for floating-point types.

#include "ext.h"

/// The call `bmm_quot(x, y)`
/// returns the quotient and remainder of `x` divided by `y`
/// in `qr` such that `qr.quot * y + qr.rem == x` and `qr.rem >= 0`.
/// If `y == 0` or `x` and `y` are infinite or `x` or `y` are not numbers,
/// the behavior is undefined.
__attribute__ ((__const__, __pure__))
inline type(bmm_quot_t, A) type(bmm_quot, A)(A const x, A const y) {
  A const q = TRUNCA(x / y);
  A const r = FMODA(x, y);
  A const s = r >= 0 ? 0 : y < 0 ? -1 : 1;
  type(bmm_quot_t, A) const qr = {.quot = q - s, .rem = r + s * y};

  return qr;
}

/// The call `bmm_abs(x)`
/// returns the absolute value of `x`.
/// If `x` is not a number, the behavior is undefined.
__attribute__ ((__const__, __pure__))
inline A type(bmm_abs, A)(A const x) {
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
inline A type(bmm_wrap, A)(A const x, A const a, A const b) {
#ifndef DEBUG
  dynamic_assert(b > a, "Invalid argument");
#endif

  return type(bmm_quot, A)(x - a, b - a).rem + a;
}

/// The call `bmm_uwrap(x, b)`
/// is equivalent to `bmm_wrap(x, 0, b)`.
/// If `b <= 0` or `x` is infinite or `x` or `b` are not numbers,
/// the behavior is undefined.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A type(bmm_uwrap, A)(A const x, A const b) {
#ifndef DEBUG
  dynamic_assert(b > 0, "Invalid argument");
#endif

  return type(bmm_quot, A)(x, b).rem;
}

/// The call `bmm_resum2(x, y)`
/// returns the reciprocal sum of `x` and `y`.
/// If `x <= 0` or `y <= 0`, `x` or `y` are infinite or
/// `x` or `y` are not numbers, the behavior is undefined.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A type(bmm_resum2, A)(A const x, A const y) {
#ifndef DEBUG
  dynamic_assert(x > 0, "Invalid argument");
  dynamic_assert(y > 0, "Invalid argument");
#endif

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
inline A type(bmm_pmean2, A)(A const x, A const y, A const e) {
  return POWA((POWA(x, e) + POWA(y, e)) / 2, 1 / e);
}

/// The call `bmm_amean2(x, y)`
/// returns the arithmetic mean of `x` and `y`.
/// It is equivalent to `bmm_pmean2(x, y, 1)`.
/// If `x` or `y` are infinite or
/// `x` or `y` are not numbers, the behavior is undefined.
__attribute__ ((__const__, __pure__))
inline A type(bmm_amean2, A)(A const x, A const y) {
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
inline A type(bmm_gmean2, A)(A const x, A const y) {
#ifndef DEBUG
  dynamic_assert(x >= 0, "Invalid argument");
  dynamic_assert(y >= 0, "Invalid argument");
#endif

  return sqrt(x * y);
}

/// The call `bmm_hmean2(x, y)`
/// returns the harmonic mean of `x` and `y`.
/// It is equivalent to `bmm_pmean2(x, y, -1)`.
/// If `x <= 0` or `y <= 0`, `x` or `y` are infinite or
/// `x` or `y` are not numbers, the behavior is undefined.
__attribute__ ((__const__, __pure__))
inline A type(bmm_hmean2, A)(A const x, A const y) {
  return 2 * type(bmm_resum2, A)(x, y);
}
