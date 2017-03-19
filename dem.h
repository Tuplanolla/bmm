// Discrete element method with some assumptions.
#ifndef BMM_DEM_H
#define BMM_DEM_H

#include "conf.h"
#include "ext.h"
#include <stdbool.h>
#include <stddef.h>

struct bmm_dem_opts {
  size_t ncell[2];
  size_t nbin;
  size_t npart;
  size_t nstep;
  // Particle size mean.
  size_t rmean;
  // Particle size standard deviation.
  size_t rstd;
};

// A terrible name for constant particle data.
struct bmm_dem_partc {
  double mass;
};

struct bmm_dem_part {
  double rrad;
  double mass; // TODO No!
  struct {
    double r[2];
    double v[2];
    double f[2];
  } lin;
  struct {
    double alpha;
    double omega;
    double tau;
  } ang;
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

struct bmm_dem_buf {
  struct bmm_dem_neigh neigh;
  struct bmm_dem_partc partcs[BMM_PART_MAX];
  struct bmm_dem_part parts[BMM_PART_MAX];
};

struct bmm_dem {
  struct bmm_dem_opts opts;
  size_t istep;
  double tstep;
  double rexts[2];
  // TODO Initial value system goes here.
  // TODO Force scheme goes here.
  void (* forcesch)(struct bmm_dem*);
  // TODO Integration scheme goes here.
  void (* intsch)(struct bmm_dem*);
  // TODO Measurement system goes here.
  // TODO Nearest neighbor system goes here.
  bool dblbuf;
  bool _pad0[7];
  union {
    struct {
      struct bmm_dem_buf bufs[2];
      struct bmm_dem_buf* active;
      struct bmm_dem_buf* passive;
    } bufs;
    struct bmm_dem_buf buf;
  } data;
};

// The call `bmm_dem_getbuf(dem)`
// returns the active read and write buffer of the simulation `dem`
// whether it is single-buffered or double-buffered.
__attribute__ ((__nonnull__))
inline struct bmm_dem_buf* bmm_dem_getbuf(struct bmm_dem* const dem) {
  return dem->dblbuf ? dem->data.bufs.active : &dem->data.buf;
}

// The call `bmm_dem_getrbuf(dem)`
// returns the active read buffer of the simulation `dem`
// whether it is single-buffered or double-buffered.
__attribute__ ((__nonnull__))
inline struct bmm_dem_buf const* bmm_dem_getrbuf(
    struct bmm_dem const* const dem) {
  return dem->dblbuf ? dem->data.bufs.active : &dem->data.buf;
}

// The call `bmm_dem_getwbuf(dem)`
// returns the passive write buffer of the simulation `dem`
// whether it is single-buffered or double-buffered.
__attribute__ ((__nonnull__))
inline struct bmm_dem_buf* bmm_dem_getwbuf(struct bmm_dem* const dem) {
  return dem->dblbuf ? dem->data.bufs.passive : &dem->data.buf;
}

// The call `bmm_dem_swapbuf(dem)`
// swaps the active and passive buffers of the simulation `dem`.
// This is necessary for operations that are not structure-preserving.
__attribute__ ((__nonnull__))
inline void bmm_dem_swapbuf(struct bmm_dem* const dem) {
  if (dem->dblbuf) {
    struct bmm_dem_buf* const buf = dem->data.bufs.passive;
    dem->data.bufs.passive = dem->data.bufs.active;
    dem->data.bufs.active = buf;
  }
}

__attribute__ ((__nonnull__))
void bmm_dem_defopts(struct bmm_dem_opts*);

__attribute__ ((__nonnull__))
void bmm_dem_defpart(struct bmm_dem_part*);

__attribute__ ((__nonnull__))
void bmm_dem_def(struct bmm_dem*, struct bmm_dem_opts const*);

__attribute__ ((__nonnull__))
bool bmm_dem_run(struct bmm_dem_opts const*);

#endif
