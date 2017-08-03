/// Half-open intervals.

#ifndef BMM_IVAL_H
#define BMM_IVAL_H

#include <stdbool.h>

#include "common.h"
#include "ext.h"
#include "fp.h"

/// The call `bmm_ival_midpoint(a)`
/// returns the midpoint of the interval `a`.
__attribute__ ((__const__, __pure__))
inline double bmm_ival_midpoint(double const *const a) {
  return type(bmm_amean2, double)(a[0], a[1]);
}

/// The call `bmm_ival_length(a)`
/// returns the length of the interval `a`.
__attribute__ ((__const__, __pure__))
inline double bmm_ival_length(double const *const a) {
  return a[1] - a[0];
}

/// The call `bmm_ival_inside(a, x)`
/// checks whether the number `x` is inside the interval `a`.
__attribute__ ((__const__, __pure__))
inline bool bmm_ival_inside(double const *const a, double const x) {
  return x >= a[0] && x < a[1];
}

#endif
