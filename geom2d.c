#include "geom2d.h"

extern inline double bmm_geom2d_norm2(double const*);

extern inline double bmm_geom2d_norm(double const*);

extern inline double bmm_geom2d_angle(double const*);

extern inline void bmm_geom2d_to_polar2(double*, double const*);

extern inline void bmm_geom2d_to_polar(double*, double const*);

extern inline void bmm_geom2d_from_polar(double*, double const*);

extern inline void bmm_geom2d_from_polar2(double*, double const*);

extern inline void bmm_geom2d_pdiff(double*,
    double const*, double const*, double const*);

extern inline double bmm_geom2d_pdist2(double const*, double const*,
    double const*);

extern inline double bmm_geom2d_pdist(double const*, double const*,
    double const*);

extern inline double bmm_geom2d_pangle(double const*, double const*,
    double const*);
