#include <stdbool.h>

#include "geom2d.h"

extern inline double bmm_geom2d_norm2(double const*);

extern inline double bmm_geom2d_norm(double const*);

extern inline double bmm_geom2d_dir(double const*);

extern inline void bmm_geom2d_add(double* restrict,
    double const* restrict, double const* restrict);

extern inline void bmm_geom2d_scale(double* restrict,
    double const* restrict, double);

extern inline void bmm_geom2d_normal(double* restrict,
    double const* restrict);

extern inline void bmm_geom2d_rperpr(double* restrict,
    double const* restrict);

extern inline void bmm_geom2d_rperp(double* restrict,
    double const* restrict);

extern inline void bmm_geom2d_lperpr(double* restrict,
    double const* restrict);

extern inline void bmm_geom2d_lperp(double* restrict,
    double const* restrict);

extern inline void bmm_geom2d_to_polar2(double* restrict,
    double const* restrict);

extern inline void bmm_geom2d_to_polar(double* restrict,
    double const* restrict);

extern inline void bmm_geom2d_from_polar(double* restrict,
    double const* restrict);

extern inline void bmm_geom2d_from_polar2(double* restrict,
    double const* restrict);

extern inline void bmm_geom2d_diff(double* restrict,
    double const* restrict, double const* restrict);

extern inline double bmm_geom2d_dist2(double const* restrict,
    double const* restrict);

extern inline double bmm_geom2d_dist(double const* restrict,
    double const* restrict);

extern inline double bmm_geom2d_angle(double const* restrict,
    double const* restrict);

extern inline void bmm_geom2d_pdiff(double* restrict,
    double const* restrict, double const* restrict, double const* restrict);

extern inline double bmm_geom2d_pdist2(double const* restrict,
    double const* restrict, double const* restrict);

extern inline double bmm_geom2d_pdist(double const* restrict,
    double const* restrict, double const* restrict);

extern inline double bmm_geom2d_pangle(double const* restrict,
    double const* restrict, double const* restrict);

extern inline void bmm_geom2d_cpdiff(double* restrict,
    double const* restrict, double const* restrict,
    bool const* restrict, double const* restrict);

extern inline double bmm_geom2d_cpdist2(double const* restrict,
    double const* restrict, bool const* restrict, double const* restrict);

extern inline double bmm_geom2d_cpdist(double const* restrict,
    double const* restrict, bool const* restrict, double const* restrict);

extern inline double bmm_geom2d_cpangle(double const* restrict,
    double const* restrict, bool const* restrict, double const* restrict);
