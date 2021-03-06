#include <stdbool.h>

#include "geom2d.h"

extern inline double bmm_geom2d_dot(double const *restrict,
    double const *restrict);

extern inline double bmm_geom2d_norm2(double const *);

extern inline double bmm_geom2d_norm(double const *);

extern inline double bmm_geom2d_dir(double const *);

extern inline double bmm_geom2d_redir(double);

extern inline void bmm_geom2d_add(double *restrict,
    double const *restrict, double const *restrict);

extern inline void bmm_geom2d_addto(double *restrict, double const *restrict);

extern inline void bmm_geom2d_scale(double *restrict,
    double const *restrict, double);

extern inline void bmm_geom2d_normal(double *restrict,
    double const *restrict);

extern inline void bmm_geom2d_project(double *restrict,
    double const *restrict, double const *restrict);

extern inline void bmm_geom2d_rperp(double *restrict,
    double const *restrict);

extern inline void bmm_geom2d_lperp(double *restrict,
    double const *restrict);

extern inline void bmm_geom2d_to_polar2(double *restrict,
    double const *restrict);

extern inline void bmm_geom2d_to_polar(double *restrict,
    double const *restrict);

extern inline void bmm_geom2d_from_polar(double *restrict,
    double const *restrict);

extern inline void bmm_geom2d_from_polar2(double *restrict,
    double const *restrict);

extern inline void bmm_geom2d_diff(double *restrict,
    double const *restrict, double const *restrict);

extern inline void bmm_geom2d_diffto(double *restrict,
    double const *restrict);

extern inline double bmm_geom2d_dist2(double const *restrict,
    double const *restrict);

extern inline double bmm_geom2d_dist(double const *restrict,
    double const *restrict);

extern inline double bmm_geom2d_angle(double const *restrict,
    double const *restrict);

extern inline void bmm_geom2d_pdiff(double *restrict,
    double const *restrict, double const *restrict, double const *restrict);

extern inline double bmm_geom2d_pdist2(double const *restrict,
    double const *restrict, double const *restrict);

extern inline double bmm_geom2d_pdist(double const *restrict,
    double const *restrict, double const *restrict);

extern inline double bmm_geom2d_pangle(double const *restrict,
    double const *restrict, double const *restrict);

extern inline void bmm_geom2d_cpdiff(double *restrict,
    double const *restrict, double const *restrict,
    double const *restrict, bool const *restrict);

extern inline double bmm_geom2d_cpdist2(double const *restrict,
    double const *restrict, double const *restrict, bool const *restrict);

extern inline double bmm_geom2d_cpdist(double const *restrict,
    double const *restrict, double const *restrict, bool const *restrict);

extern inline double bmm_geom2d_cpangle(double const *restrict,
    double const *restrict, double const *restrict, bool const *restrict);

extern inline void bmm_geom2d_refl(double *restrict,
    double const *restrict, double const *restrict, int);

extern inline void bmm_geom2d_reflto(double *restrict,
    double const *restrict, int);

extern inline double bmm_geom2d_shellvol(double const *restrict,
    double, double const *restrict, bool const *);
