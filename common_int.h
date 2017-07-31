/// Common operations for integer types.

#include <stdbool.h>

#include "ext.h"

/// The call `bmm_even(x)`
/// checks whether `x` is even.
__attribute__ ((__const__, __pure__))
inline bool type(bmm_even, A)(A const x) {
  return x % 2 == 0;
}

/// The call `bmm_odd(x)`
/// checks whether `x` is odd.
__attribute__ ((__const__, __pure__))
inline bool type(bmm_odd, A)(A const x) {
  return x % 2 != 0;
}

/// The call `bmm_tamean2(x, y)`
/// returns the truncated arithmetic mean of `x` and `y`.
/// Overflows are impossible both internally and externally.
__attribute__ ((__const__, __pure__))
inline A type(bmm_tamean2, A)(A const x, A const y) {
  // TODO Not like this.

  A const z = (x / 2 + y / 2);

  if ((x % 2 == 0) && (y % 2 == 0))
    return z;

  if ((x % 2 == 1) && (y % 2 == 1))
    return z + (z >= 0); // z + 1

  if ((x % 2 == -1) && (y % 2 == -1))
    return z - (z <= 0); // z - 1

  if ((x % 2 == 0 && y % 2 == 1) ||
      (x % 2 == 1 && y % 2 == 0))
    return z + (z < 0);

  if ((x % 2 == 1 && y % 2 == -1) ||
      (x % 2 == -1 && y % 2 == 1))
    return z;

  if ((x % 2 == -1 && y % 2 == 0) ||
      (x % 2 == 0 && y % 2 == -1))
    return z - (z > 0);

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
