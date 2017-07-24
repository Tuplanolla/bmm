/// Common operations for floating-point types.

#include "ext.h"

/// The call `bmm_quot(x, y)`
/// returns the quotient and remainder of `x` divided by `y`
/// in `qr` such that `qr.quot * y + qr.rem == x` and `qr.rem >= 0`.
/// If `y == 0` or `x` and `y` are infinite or `x` or `y` are not numbers,
/// the behavior is undefined.
__attribute__ ((__const__, __pure__))
inline type(bmm_quot_t, A) type(bmm_quot, A)(A const x, A const y) {
  A const q = trunc(x / y);
  A const r = fmod(x, y);
  A const s = r >= 0 ? 0 : y < 0 ? -1 : 1;
  type(bmm_quot_t, A) const qr = {.quot = q - s, .rem = r + s * y};

  return qr;
}

// TODO Wrapping can be extended for signed and floating-point `A`.

/// The call `bmm_wrap(x, a, b)`
/// finds such `y` that `a <= y < b`
/// by shifting `x` by the appropriate number of multiples of `b - a`.
/// If `b <= a`, the behavior is undefined.
/// Overflows are impossible both internally and externally.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A type(bmm_wrap, A)(A const x, A const a, A const b) {
#ifndef DEBUG
  dynamic_assert(b > a, "Invalid argument");
#endif

  A const c = b - a;

  A y = x;

  if (y < a)
    do
      y += c;
    while (y < a);
  else if (y >= b)
    do
      y -= c;
    while (y >= b);

  return y;
}
