#ifndef BMM_DEM_H
#define BMM_DEM_H

#include <stdbool.h>
#include <stddef.h>

#define DIM_MAX ((size_t) 2)
#define BIN_MAX ((size_t) 1024)
#define PART_MAX ((size_t) 16384)
#define STEP_MAX ((size_t) 16777216)

struct bmm_opts {
  size_t ndim;
  size_t nbin;
  size_t npart;
  size_t nstep;
};

struct bmm_part {
  double rpos[DIM_MAX];
};

struct bmm_state {
  struct bmm_opts opts;
  size_t istep;
  double rexts[DIM_MAX];
  struct bmm_part parts[PART_MAX];
};

bool bmm_rundem(struct bmm_opts const*);

#endif
