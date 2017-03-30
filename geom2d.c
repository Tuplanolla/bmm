#include "geom2d.h"

extern inline double bmm_geom2d_norm2(double const*);

extern inline double bmm_geom2d_norm(double const*);

extern inline double bmm_geom2d_angle(double const*);

extern inline void bmm_geom2d_add(double*, double const*, double const*);

extern inline void bmm_geom2d_scale(double*, double const*, double);

extern inline void bmm_geom2d_normal(double*, double const*);

extern inline void bmm_geom2d_rperpr(double* restrict, double const* restrict);

extern inline void bmm_geom2d_rperp(double*, double const*);

extern inline void bmm_geom2d_lperpr(double* restrict, double const* restrict);

extern inline void bmm_geom2d_lperp(double*, double const*);

extern inline void bmm_geom2d_to_polar2(double*, double const*);

extern inline void bmm_geom2d_to_polar(double*, double const*);

extern inline void bmm_geom2d_from_polar(double*, double const*);

extern inline void bmm_geom2d_from_polar2(double*, double const*);

// TODO Organize these.

extern inline void bmm_geom2d_diff(double*, double const*, double const*);

extern inline double bmm_geom2d_dist2(double const*, double const*);

extern inline double bmm_geom2d_dist(double const*, double const*);

extern inline void bmm_geom2d_pdiff(double*,
    double const*, double const*, double const*);

extern inline double bmm_geom2d_pdist2(double const*, double const*,
    double const*);

extern inline double bmm_geom2d_pdist(double const*, double const*,
    double const*);

extern inline double bmm_geom2d_pangle(double const*, double const*,
    double const*);
