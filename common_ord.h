/// Common operations for ordered types.

#include <stdbool.h>
#include <stddef.h>

#include "cpp.h"
#include "ext.h"

/// The call `bmm_cmp(x, y)`
/// returns the comparison of `x` and `y`, which is
///
/// * `-1` if `x < y`,
/// * `1` if `x > y` and
/// * `0` otherwise.
///
/// This is useful with `bmm_hsort` for example.
__attribute__ ((__const__, __pure__))
inline int type(bmm_cmp, A)(A const x, A const y) {
  return x < y ? -1 : x > y ? 1 : 0;
}

/// The call `bmm_clamp(x, a, b)`
/// finds such `y` that `a <= y <= b`
/// by shifting `x` by the smallest possible amount.
/// If `b < a` or `x` is infinite or `x`, `a` or `b` are not numbers,
/// the behavior is undefined.
/// Overflows are impossible both internally and externally.
__attribute__ ((__const__, __pure__))
inline A type(bmm_clamp, A)(A const x, A const a, A const b) {
  dynamic_assert(b >= a, "Invalid argument");

  return x < a ? a : x > b ? b : x;
}

/// The call `bmm_min(x, y)`
/// returns the lesser of `x` and `y`.
/// If `x` or `y` are not numbers, the behavior is undefined.
__attribute__ ((__const__, __pure__))
inline A type(bmm_min, A)(A const x, A const y) {
  return x < y ? x : y;
}

/// The call `bmm_max(x, y)`
/// returns the lesser of `x` and `y`.
/// If `x` or `y` are not numbers, the behavior is undefined.
__attribute__ ((__const__, __pure__))
inline A type(bmm_max, A)(A const x, A const y) {
  return x > y ? x : y;
}
