/// Common operations for numeric types.

#include <stddef.h>

#include "ext.h"

/// This structure holds the quotient and remainder of division
/// in unspecified order.
typedef struct {
  A quot;
  A rem;
} type(bmm_quot_t, A);

/// The call `bmm_sgn(x)`
/// returns the sign of `x`.
/// Overflows are impossible both internally and externally.
__attribute__ ((__const__, __pure__))
inline int type(bmm_sgn, A)(A const x) {
  return type(bmm_cmp, A)(x, 0);
}

/// The call `bmm_abs(x)`
/// returns the absolute value of `x`.
/// Overflows are impossible internally but possible externally.
__attribute__ ((__const__, __pure__))
inline A type(bmm_abs, A)(A const x) {
  return x < 0 ? -x : x;
}

/// The call `bmm_pow(x, e)`
/// returns `x` raised to the power of `e`.
/// If `x == 0` and `e == 0`, the result is one.
/// Overflows are impossible internally but possible externally.
__attribute__ ((__const__, __pure__))
inline A type(bmm_pow, A)(A const x, size_t const e) {
  A y = 1;

  A m = x;
  for (size_t p = e; p > 0; m *= m, p >>= 1)
    if ((p & 1) != 0)
      y *= m;

  return y;

  // The following implementation is less complicated,
  // but slower for large powers.
  // A y = 1;
  //
  // for (size_t i = 0; i < e; ++i)
  //   y *= x;
  //
  // return y;
}
