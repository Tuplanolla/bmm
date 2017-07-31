/// Common operations for signed integer types.

#include "ext.h"

/// The call `bmm_quot(x, y)`
/// returns the quotient and remainder of `x` divided by `y`
/// in `qr` such that `qr.quot * y + qr.rem == x` and `qr.rem >= 0`.
/// If `y == 0`, the behavior is undefined.
/// Overflows are impossible both internally and externally.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline type(bmm_quot_t, A) type(bmm_quot, A)(A const x, A const y) {
#ifndef DEBUG
  dynamic_assert(y != 0, "Invalid argument");
#endif

  A const q = x / y;
  A const r = x % y;
  A const s = r >= 0 ? 0 : y < 0 ? -1 : 1;
  type(bmm_quot_t, A) const qr = {.quot = q - s, .rem = r + s * y};

  return qr;
}

/// The call `bmm_abs(x)`
/// returns the absolute value of `x`.
/// If `-x` is not representable, the behavior is undefined.
/// Overflows are impossible internally but possible externally.
__attribute__ ((__const__, __pure__))
inline A type(bmm_abs, A)(A const x) {
  return x < 0 ? -x : x;
}

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

  // The condition below guarantees
  // that the body of the loop at the end is executed at most twice.

  A const m = type(bmm_abs, A)(x / 2);
  if (a >= 0 || b <= 0 || (a > -m && b < m)) {
    A const c = b - a;
    A const r = type(bmm_quot, A)(x, c).rem;
    A const s = type(bmm_quot, A)(a, c).rem;

    return (r >= s ? r - s : c - (s - r)) + a;
  }

  A y = x;

  if (y < a)
    do {
      y += b;
      y -= a;
    } while (y < a);
  else if (y >= b)
    do {
      y -= b;
      y += a;
    } while (y >= b);

  return y;
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

/// The call `bmm_tamean2(x, y)`
/// returns the truncated arithmetic mean of `x` and `y`.
/// Overflows are impossible both internally and externally.
__attribute__ ((__const__, __pure__))
inline A type(bmm_tamean2, A)(A const x, A const y) {
  A const z = x / 2 + y / 2;

  // This function was derived with the help
  // of the following Espresso specification.
  //
  //     .i 4
  //     .o 6
  //     .ilb sx mx sy my
  //     .ob sl ml se me sg mg
  //     01 01 -0 01 01
  //     01 -0 01 -0 -0
  //     01 11 -0 -0 -0
  //     -0 01 01 -0 -0
  //     -0 -0 -0 -0 -0
  //     -0 11 -0 -0 11
  //     11 01 -0 -0 -0
  //     11 -0 -0 -0 11
  //     11 11 11 11 -0
  //     .e
  //
  // The input variables are defined as
  //
  // * `sx = x % 2 < 0`,
  // * `mx = x % 2 != 0`,
  // * `sy = y % 2 < 0` and
  // * `my = y % 2 != 0`.
  //
  // The output variables are then used to define
  //
  // * `l = (sl ? -1 : 1) * (ml ? 1 : 0) * (z < 0 ? 1 : 0)`,
  // * `e = (se ? -1 : 1) * (me ? 1 : 0) * (z == 0 ? 1 : 0)` and
  // * `g = (sg ? -1 : 1) * (mg ? 1 : 0) * (z > 0 ? 1 : 0)`.
  //
  // They form the result `z + l + e + g`.

  switch (x % 2 + y % 2) {
    case -2:
      return z - (z < 0 ? 1 : 0) - (z == 0 ? 1 : 0);
    case -1:
      return z - (z > 0 ? 1 : 0);
    case 0:
      return z;
    case 1:
      return z + (z < 0 ? 1 : 0);
    case 2:
      return z + (z > 0 ? 1 : 0) + (z == 0 ? 1 : 0);
  }

  return z;

  // The following implementation is less complicated,
  // but susceptible to overflowing.
  // return (x + y) / 2;
}

/// The call `bmm_famean2(x, y)`
/// returns the floored arithmetic mean of `x` and `y`.
/// Overflows are impossible both internally and externally.
__attribute__ ((__const__, __pure__))
inline A type(bmm_famean2, A)(A const x, A const y) {
  type(bmm_quot_t, A) const qrx = type(bmm_quot, A)(x, 2);
  type(bmm_quot_t, A) const qry = type(bmm_quot, A)(y, 2);

  return (qrx.quot + qry.quot) + qrx.rem * qry.rem;

  // The following implementation is less complicated,
  // but susceptible to overflowing.
  // return type(bmm_quot, A)(x + y, 2).quot;
}
