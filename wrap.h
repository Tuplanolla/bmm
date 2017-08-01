/// Wrappers for standard library operations.

#ifndef BMM_WRAP_H
#define BMM_WRAP_H

#include <math.h>

#include "alias.h"
#include "ext.h"

/// The call `bmm_trunc(x)`
/// returns the truncated value of `x`.
/// If `x` is infinite or `x` is not a number, the behavior is undefined.
__attribute__ ((__const__, __pure__))
inline float type(bmm_trunc, float)(float const x) {
  return truncf(x);
}

/// The call `bmm_trunc(x)`
/// returns the truncated value of `x`.
/// If `x` is infinite or `x` is not a number, the behavior is undefined.
__attribute__ ((__const__, __pure__))
inline double type(bmm_trunc, double)(double const x) {
  return trunc(x);
}

/// The call `bmm_trunc(x)`
/// returns the truncated value of `x`.
/// If `x` is infinite or `x` is not a number, the behavior is undefined.
__attribute__ ((__const__, __pure__))
inline long double type(bmm_trunc, long_double)(long double const x) {
  return truncl(x);
}

/// The call `bmm_fmod(x, y)`
/// returns the remainder of `x` divided by `y`.
/// If `y == 0` or `x` and `y` are infinite or `x` or `y` are not numbers,
/// the behavior is undefined.
__attribute__ ((__const__, __pure__))
inline float type(bmm_fmod, float)(float const x, float const y) {
  return fmodf(x, y);
}

/// The call `bmm_fmod(x, y)`
/// returns the remainder of `x` divided by `y`.
/// If `y == 0` or `x` and `y` are infinite or `x` or `y` are not numbers,
/// the behavior is undefined.
__attribute__ ((__const__, __pure__))
inline double type(bmm_fmod, double)(double const x, double const y) {
  return fmod(x, y);
}

/// The call `bmm_fmod(x, y)`
/// returns the remainder of `x` divided by `y`.
/// If `y == 0` or `x` and `y` are infinite or `x` or `y` are not numbers,
/// the behavior is undefined.
__attribute__ ((__const__, __pure__))
inline long double type(bmm_fmod, long_double)(long double const x,
    long double const y) {
  return fmodl(x, y);
}

/// The call `bmm_pow(x, y)`
/// returns `x` raised to the power of `y`.
/// If `x == 0` and `y == 0`, the result is one.
/// If `x` or `y` are not numbers, the behavior is undefined.
__attribute__ ((__const__, __pure__))
inline float type(bmm_pow, float)(float const x, float const y) {
  return powf(x, y);
}

/// The call `bmm_pow(x, y)`
/// returns `x` raised to the power of `y`.
/// If `x == 0` and `y == 0`, the result is one.
/// If `x` or `y` are not numbers, the behavior is undefined.
__attribute__ ((__const__, __pure__))
inline double type(bmm_pow, double)(double const x, double const y) {
  return pow(x, y);
}

/// The call `bmm_pow(x, y)`
/// returns `x` raised to the power of `y`.
/// If `x == 0` and `y == 0`, the result is one.
/// If `x` or `y` are not numbers, the behavior is undefined.
__attribute__ ((__const__, __pure__))
inline long double type(bmm_pow, long_double)(long double const x,
    long double const y) {
  return powl(x, y);
}

#endif
