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
