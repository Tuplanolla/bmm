// Discrete element method with some assumptions.
#ifndef BMM_DEM_H
#define BMM_DEM_H

#include "conf.h"
#include "ext.h"
#include <stdbool.h>
#include <stddef.h>

struct bmm_dem_opts {
  size_t ndim;
  size_t nbin;
  size_t npart;
  size_t nstep;
};

struct bmm_dem_part {
  double rrad;
  double arot;
  double rpos[BMM_DIM_MAX];
};

struct bmm_dem {
  struct bmm_dem_opts opts;
  size_t istep;
  double rexts[BMM_DIM_MAX];
  struct bmm_dem_part parts[BMM_PART_MAX];
};

__attribute__ ((__nonnull__))
void bmm_dem_defopts(struct bmm_dem_opts*);

__attribute__ ((__nonnull__))
void bmm_dem_defpart(struct bmm_dem_part*);

__attribute__ ((__nonnull__))
void bmm_dem_def(struct bmm_dem*, struct bmm_dem_opts const*);

__attribute__ ((__nonnull__))
bool bmm_dem_run(struct bmm_dem_opts const*);

#endif
