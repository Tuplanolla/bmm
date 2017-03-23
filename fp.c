#include "fp.h"
#include <stddef.h>

extern inline bmm_fp_div_t bmm_fp_div(double, double);

extern inline int bmm_fp_cmp(double, double);

extern inline double bmm_fp_identity(double);

extern inline double bmm_fp_constant(double, double);

extern inline double bmm_fp_zero(double);

extern inline double bmm_fp_one(double);

extern inline double bmm_fp_midpoint(double, double);

extern inline double bmm_fp_sq(double);

extern inline double bmm_fp_rt(double, double);

extern inline double bmm_fp_log(double, double);

extern inline double bmm_fp_clamp(double, double, double);

extern inline double bmm_fp_sclamp(double, double);

extern inline double bmm_fp_uclamp(double, double);

extern inline double bmm_fp_wrap(double, double, double);

extern inline double bmm_fp_swrap(double, double);

extern inline double bmm_fp_uwrap(double, double);

extern inline double bmm_fp_lerp(double, double, double, double, double);

extern inline double bmm_fp_lorp(double, double, double, double, double);

extern inline double bmm_fp_ballvol(double, size_t);

extern inline double bmm_fp_ballsurf(double, size_t);

extern inline double bmm_fp_norm2(double const*);

extern inline double bmm_fp_norm(double const*);

extern inline double bmm_fp_angle(double const*);

extern inline void bmm_fp_to_polar2(double*, double const*);

extern inline void bmm_fp_to_polar(double*, double const*);

extern inline void bmm_fp_from_polar(double*, double const*);

extern inline void bmm_fp_from_polar2(double*, double const*);

extern inline void bmm_fp_pdiff(double*,
    double const*, double const*, double const*);

extern inline double bmm_fp_pdist2(double const*, double const*,
    double const*);

extern inline double bmm_fp_pdist(double const*, double const*,
    double const*);

extern inline double bmm_fp_pangle(double const*, double const*,
    double const*);
