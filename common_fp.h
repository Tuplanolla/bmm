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
/// If `b <= 0`, the behavior is undefined.
/// Overflows are impossible both internally and externally.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A type(bmm_uwrap, A)(A const x, A const b) {
#ifndef DEBUG
  dynamic_assert(b > 0, "Invalid argument");
#endif

  return type(bmm_quot, A)(x, b).rem;
}
