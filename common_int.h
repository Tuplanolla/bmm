/// Common operations for integer types.

#include <stddef.h>

#include "ext.h"

// TODO Six different division modes seem to exist, because people are morons.

/// The call `bmm_sgn(x)`
/// returns the sign of `x`.
__attribute__ ((__const__, __pure__))
inline int type(bmm_sgn, A)(A const x) {
  return type(bmm_cmp, A)(x, 0);
}

/// This structure holds the quotient and remainder of division
/// in unspecified order.
typedef struct {
  A quot;
  A rem;
} type(bmm_quot_t, A);

/// The call `bmm_quot(x, y)`
/// returns the quotient and remainder of `x` divided by `y`
/// in `qr` such that `y * qr.quot + qr.rem == x` and
/// `bmm_sgn(qr.rem) * bmm_sgn(x) != -1`.
/// Overflows are impossible both internally and externally.
__attribute__ ((__const__, __pure__))
inline type(bmm_quot_t, A) type(bmm_quot, A)(A const x, A const y) {
  type(bmm_quot_t, A) const qr = {.quot = x / y, .rem = x % y};

  return qr;
}

/// This structure holds the division and modulo of division
/// in unspecified order.
typedef struct {
  A div;
  A mod;
} type(bmm_div_t, A);

/// The call `bmm_div(x, y)`
/// returns the division and modulo of `x` divided by `y`
/// in `dm` such that `y * dm.div + dm.mod == x` and
/// `bmm_sgn(dm.rem) * bmm_sgn(y) != -1`.
/// Overflows are impossible both internally and externally.
__attribute__ ((__const__, __pure__))
inline type(bmm_div_t, A) type(bmm_div, A)(A const x, A const y) {
  if (y < 0)
    return type(bmm_div, A)(x, -y);

  type(bmm_quot_t, A) const qr = type(bmm_quot, A)(x, y);

  type(bmm_div_t, A) dm = {.div = qr.quot, .mod = qr.rem};
  if (qr.rem < 0) {
    dm.div -= 1;
    dm.mod += y;
  }

  return dm;
}
