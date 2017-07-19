#include "ext.h"

/// The call `bmm_swap(x, y)`
/// exchanges `x` with `y`.
__attribute__ ((__nonnull__))
inline void inst(bmm_swap, T)(T *restrict const x, T *restrict const y) {
  T const tmp = *x;
  *x = *y;
  *y = tmp;
}
