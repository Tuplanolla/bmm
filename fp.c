#include <stddef.h>

#include "fp.h"

extern inline bmm_fp_div_t bmm_fp_div(double, double);

extern inline int bmm_fp_cmp(double, double);

extern inline int bmm_fp_lexcmp(double const* restrict,
    double const* restrict, size_t);

extern inline double bmm_fp_identity(double);

extern inline double bmm_fp_constant(double, double);

extern inline double bmm_fp_zero(double);

extern inline double bmm_fp_one(double);

extern inline double bmm_fp_midpoint(double, double);

extern inline double bmm_fp_sq(double);

extern inline double bmm_fp_cb(double);

extern inline double bmm_fp_rt(double, size_t);

extern inline double bmm_fp_log(double, double);

extern inline double bmm_fp_clamp(double, double, double);

extern inline double bmm_fp_sclamp(double, double);

extern inline double bmm_fp_uclamp(double, double);

extern inline double bmm_fp_wrap(double, double, double);

extern inline double bmm_fp_swrap(double, double);

extern inline double bmm_fp_uwrap(double, double);

extern inline double bmm_fp_sum(double const*, size_t);

extern inline double bmm_fp_prod(double const*, size_t);

extern inline double bmm_fp_min(double const*, size_t);

extern inline double bmm_fp_max(double const*, size_t);

extern inline double bmm_fp_lfold(double (*)(double, double, void*),
    double const* restrict, size_t, double, void* restrict);

extern inline double bmm_fp_rfold(double (*)(double, double, void*),
    double const* restrict, size_t, double, void* restrict);

extern inline double bmm_fp_lerp(double, double, double, double, double);

extern inline double bmm_fp_lorp(double, double, double, double, double);

extern inline double bmm_fp_percent(double, double);

extern inline double bmm_fp_pow(double, size_t);

extern inline double bmm_fp_fact(size_t, size_t);

extern inline size_t bmm_fp_iclerp(double, double, double, size_t, size_t);
