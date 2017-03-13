// Discrete element method with some assumptions.
#ifndef BMM_DEM_H
#define BMM_DEM_H

#include "ext.h"
#include <stdbool.h>
#include <stddef.h>

#define DIM_MAX ((size_t) 2)
#define BIN_MAX ((size_t) 1024)
// #define PART_MAX ((size_t) 16384)
#define PART_MAX ((size_t) 8)
#define STEP_MAX ((size_t) 16777216)

struct bmm_opts {
  size_t ndim;
  size_t nbin;
  size_t npart;
  size_t nstep;
};

struct bmm_part {
  double rrad;
  double rpos[DIM_MAX];
};

struct bmm_state {
  struct bmm_opts opts;
  size_t istep;
  double rexts[DIM_MAX];
  struct bmm_part parts[PART_MAX];
};

__attribute__ ((__nonnull__))
void bmm_defopts(struct bmm_opts*);

__attribute__ ((__nonnull__))
void bmm_defpart(struct bmm_part*);

__attribute__ ((__nonnull__))
void bmm_defstate(struct bmm_state*);

__attribute__ ((__nonnull__))
bool bmm_rundem(struct bmm_opts const*);

#endif
