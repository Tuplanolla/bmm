#include <stddef.h>

#include "fp.h"

extern inline int bmm_fp_lexcmp(double const *restrict,
    double const *restrict, size_t);

extern inline double bmm_fp_rt(double, size_t);

extern inline double bmm_fp_log(double, double);

extern inline double bmm_fp_clamp(double, double, double);

extern inline double bmm_fp_sclamp(double, double);

extern inline double bmm_fp_uclamp(double, double);

extern inline double bmm_fp_wrap(double, double, double);

extern inline double bmm_fp_swrap(double, double);

extern inline double bmm_fp_uwrap(double, double);

extern inline double bmm_fp_lerp(double, double, double, double, double);

extern inline double bmm_fp_lorp(double, double, double, double, double);

extern inline double bmm_fp_percent(double, double);

extern inline double bmm_fp_min(double const *, size_t);

extern inline double bmm_fp_max(double const *, size_t);

extern inline size_t bmm_fp_iclerp(double, double, double, size_t, size_t);
