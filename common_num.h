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
__attribute__ ((__const__, __pure__))
inline int type(bmm_sgn, A)(A const x) {
  return type(bmm_cmp, A)(x, 0);
}
