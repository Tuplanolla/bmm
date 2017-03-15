// Discrete element method with some assumptions.
#ifndef BMM_DEM_H
#define BMM_DEM_H

#include "ext.h"
#include <stdbool.h>
#include <stddef.h>

#define BMM_DIM_MAX ((size_t) 2)
#define BMM_BIN_MAX ((size_t) 1024)
// #define BMM_PART_MAX ((size_t) 16384)
#define BMM_PART_MAX ((size_t) 8)
#define BMM_STEP_MAX ((size_t) 16777216)

struct bmm_dem_opts {
  size_t ndim;
  size_t nbin;
  size_t npart;
  size_t nstep;
};

struct bmm_part {
  double rrad;
  double arot;
  double rpos[BMM_DIM_MAX];
};

struct bmm_state {
  struct bmm_dem_opts opts;
  size_t istep;
  double rexts[BMM_DIM_MAX];
  struct bmm_part parts[BMM_PART_MAX];
};

__attribute__ ((__nonnull__))
void bmm_defopts(struct bmm_dem_opts*);

__attribute__ ((__nonnull__))
void bmm_defpart(struct bmm_part*);

__attribute__ ((__nonnull__))
void bmm_defstate(struct bmm_state*);

__attribute__ ((__nonnull__))
bool bmm_dem_run(struct bmm_dem_opts const*);

#endif
