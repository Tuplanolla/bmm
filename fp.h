// Floating-point operations.
#ifndef BMM_FP_H
#define BMM_FP_H

#ifndef M_PI
#define M_PI 3.141592653589793
#endif

#ifndef M_2PI
#define M_2PI 6.283185307179586
#endif

__attribute__ ((__const__, __pure__))
inline double bmm_fp_lerp(double const x,
    double const x0, double const x1, double const y0, double const y1) {
  return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
}

#endif
