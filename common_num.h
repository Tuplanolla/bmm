/// Common operations for numeric types.

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
