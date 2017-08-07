/// Wrappers for standard library operations.

#ifndef BMM_WRAP_H
#define BMM_WRAP_H

#include <math.h>

#include "alias.h"
#include "ext.h"

__attribute__ ((__const__, __pure__))
inline double type(fabs, double)(double const x) {
  return fabs(x);
}

__attribute__ ((__const__, __pure__))
inline double type(trunc, double)(double const x) {
  return trunc(x);
}

__attribute__ ((__const__, __pure__))
inline double type(fmod, double)(double const x, double const y) {
  return fmod(x, y);
}

__attribute__ ((__const__, __pure__))
inline double type(pow, double)(double const x, double const y) {
  return pow(x, y);
}

#endif
