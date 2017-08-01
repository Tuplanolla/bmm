/// Wrappers for standard library operations.

#ifndef BMM_WRAP_H
#define BMM_WRAP_H

#include <math.h>

#include "alias.h"
#include "ext.h"

__attribute__ ((__const__, __pure__))
inline float type(fabs, float)(float const x) {
  return fabsf(x);
}

__attribute__ ((__const__, __pure__))
inline double type(fabs, double)(double const x) {
  return fabs(x);
}

__attribute__ ((__const__, __pure__))
inline long double type(fabs, long_double)(long double const x) {
  return fabsl(x);
}

__attribute__ ((__const__, __pure__))
inline float type(trunc, float)(float const x) {
  return truncf(x);
}

__attribute__ ((__const__, __pure__))
inline double type(trunc, double)(double const x) {
  return trunc(x);
}

__attribute__ ((__const__, __pure__))
inline long double type(trunc, long_double)(long double const x) {
  return truncl(x);
}

__attribute__ ((__const__, __pure__))
inline float type(fmod, float)(float const x, float const y) {
  return fmodf(x, y);
}

__attribute__ ((__const__, __pure__))
inline double type(fmod, double)(double const x, double const y) {
  return fmod(x, y);
}

__attribute__ ((__const__, __pure__))
inline long double type(fmod, long_double)(long double const x,
    long double const y) {
  return fmodl(x, y);
}

__attribute__ ((__const__, __pure__))
inline float type(pow, float)(float const x, float const y) {
  return powf(x, y);
}

__attribute__ ((__const__, __pure__))
inline double type(pow, double)(double const x, double const y) {
  return pow(x, y);
}

__attribute__ ((__const__, __pure__))
inline long double type(pow, long_double)(long double const x,
    long double const y) {
  return powl(x, y);
}

#endif
