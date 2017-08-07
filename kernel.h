#ifndef BMM_KERNEL_H
/// Statistical kernels.
#define BMM_KERNEL_H

#include <math.h>

#include "common.h"
#include "ext.h"
#include "fp.h"

/// The call `bmm_kernel_rect(x)`
/// returns the value of the uniform rectangular kernel at `x`.
/// The support of this kernel is the signed unit interval.
__attribute__ ((__const__, __pure__))
inline double bmm_kernel_rect(double const x) {
  return fabs(x) < 1.0 ? 1.0 / 2.0 : 0.0;
}

/// The call `bmm_kernel_tri(x)`
/// returns the value of the triangular kernel at `x`.
/// The support of this kernel is the signed unit interval.
__attribute__ ((__const__, __pure__))
inline double bmm_kernel_tri(double const x) {
  double const absx = fabs(x);

  return absx < 1.0 ? 1.0 - absx : 0.0;
}

/// The call `bmm_kernel_epan(x)`
/// returns the value of the parabolic Epanechnikov kernel at `x`.
/// The support of this kernel is the signed unit interval.
__attribute__ ((__const__, __pure__))
inline double bmm_kernel_epan(double const x) {
  double const x2 = $(bmm_power, double)(x, 2);

  return x2 < 1.0 ? (3.0 / 4.0) * (1.0 - x2) : 0.0;
}

/// The call `bmm_kernel_biweight(x)`
/// returns the value of the quartic biweight kernel at `x`.
/// The support of this kernel is the signed unit interval.
__attribute__ ((__const__, __pure__))
inline double bmm_kernel_biweight(double const x) {
  double const x2 = $(bmm_power, double)(x, 2);

  return x2 < 1.0 ? (15.0 / 16.0) * $(bmm_power, double)(1.0 - x2, 2) : 0.0;
}

/// The call `bmm_kernel_cos(x)`
/// returns the value of the cosine kernel at `x`.
/// The support of this kernel is the signed unit interval.
__attribute__ ((__const__, __pure__))
inline double bmm_kernel_cos(double const x) {
  double const absx = fabs(x);

  return absx < 1.0 ? M_PI_4 * cos(M_PI_2 * x) : 0.0;
}

/// The call `bmm_kernel_gaussian(x)`
/// returns the value of the Gaussian kernel at `x`.
/// The support of this kernel is the entire real line.
__attribute__ ((__const__, __pure__))
inline double bmm_kernel_gaussian(double const x) {
  return (1.0 / sqrt(M_2PI)) *
    exp(-(1.0 / 2.0) * $(bmm_power, double)(x, 2));
}

/// The call `bmm_kernel_logistic(x)`
/// returns the value of the logistic kernel at `x`.
/// The support of this kernel is the entire real line.
__attribute__ ((__const__, __pure__))
inline double bmm_kernel_logistic(double const x) {
  return (1.0 / 4.0) / $(bmm_power, double)(cosh((1.0 / 2.0) * x), 2);
}

enum bmm_kernel {
  BMM_KERNEL_RECT,
  BMM_KERNEL_TRI,
  BMM_KERNEL_EPAN,
  BMM_KERNEL_BIWEIGHT,
  BMM_KERNEL_COS,
  BMM_KERNEL_GAUSSIAN,
  BMM_KERNEL_LOGISTIC
};

/// The call `bmm_kernel(k)`
/// returns the function for the kernel `k`.
__attribute__ ((__const__, __pure__))
inline double (*bmm_kernel(enum bmm_kernel const k))(double) {
  static double (*const bmm_kernels[])(double) = {
    bmm_kernel_rect,
    bmm_kernel_tri,
    bmm_kernel_epan,
    bmm_kernel_biweight,
    bmm_kernel_cos,
    bmm_kernel_gaussian,
    bmm_kernel_logistic
  };

  return bmm_kernels[k];
}

#endif
