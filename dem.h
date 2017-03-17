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
  // Particle size mean.
  size_t rmean;
  // Particle size standard deviation.
  size_t rstd;
};

struct bmm_dem_part {
  double rrad;
  double arot;
  double rpos[2];
};

struct bmm_dem_list {
  // Size prefix.
  size_t n;
  size_t i[BMM_GROUP_MAX];
};

struct bmm_dem_neigh {
  // Maximum distance for qualifying as a neighbor.
  double rmax;
  // Cell extents by dimension, expressed in divs.
  size_t ncell[2];
  // Which particles each cell contains, roughly.
  struct bmm_dem_list parts[BMM_CELL_MAX * BMM_CELL_MAX];
  // Neighbors for each particle.
  struct bmm_dem_list neighs[BMM_PART_MAX];
};

struct bmm_dem {
  struct bmm_dem_opts opts;
  size_t istep;
  double rexts[2];
  // TODO Initial value system goes here.
  // TODO Integration scheme goes here.
  void (* intsch)(struct bmm_dem*);
  // TODO Force scheme goes here.
  // TODO Measurement system goes here.
  // TODO Nearest neighbor system goes here.
  struct bmm_dem_neigh neigh;
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
