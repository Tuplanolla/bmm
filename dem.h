#ifndef BMM_DEM_H
/// Discrete element method with some assumptions.
#define BMM_DEM_H

#include <gsl/gsl_rng.h>
#include <stdbool.h>
#include <stddef.h>

#include "conf.h"
#include "ext.h"

enum bmm_dem_mode {
  BMM_DEM_BEGIN,
  BMM_DEM_SEDIMENT,
  BMM_DEM_LINK,
  BMM_DEM_BREAK,
  BMM_DEM_ACCEL,
  BMM_DEM_CRUNCH
};

struct bmm_dem_est {
  double ekinetic;
  double pvector;
  double pscalar;
};

struct bmm_dem_time {
  double tinit;
  double tsim;
};

struct bmm_dem_opts {
  // Cell extents by dimension, expressed in divs.
  size_t ncell[2];
  size_t nbin;
  // __attribute__ ((__deprecated__))
  size_t nstep;
  // Maximum distance for qualifying as a neighbor.
  double rmax;
  // TODO Total simulation time.
  struct bmm_dem_time tend;
  struct bmm_dem_time tadv;
  // Simulated time step.
  // __attribute__ ((__deprecated__))
  double tstep;
  // __attribute__ ((__deprecated__))
  double tstepcomm;
  double tcomm;
  // Link length creation multiplier.
  double linkslurp;
  // Link strength (according to Hooke's law).
  double klink;
  // Drift leeway velocity.
  double vleeway;
  // Cohesive force multiplier during sedimentation.
  double fcohes;
  // Accelerative force during crunching.
  double faccel;
  // Acceleration for sedimentation (gravity).
  double gravy[2];
  // Young's modulus $Y$.
  double ymodul;
  // Elasticity $\gamma$.
  double yelast;
  // Nonphysical ambient damping factor.
  double damp;
  // Particle size mean.
  double rmean;
  // Particle size standard deviation.
  double rsd;
};

// A terrible name for constant particle data.
struct bmm_dem_partc {
  double rrad;
  double mass;
  double moi;
  // Fixed or moving.
  bool free;
};

struct bmm_dem_part {
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

struct bmm_dem_thingy {
  // Neighbor index.
  size_t i;
  // Compression history.
  double x[2];
};

struct bmm_dem_listy {
  size_t n;
  struct bmm_dem_thingy thingy[BMM_GROUP_MAX];
};

struct bmm_dem_neigh {
  // Next scheduled update.
  double tnext;
  // Neighbors for each particle
  // plus room for carrying garbage.
  struct bmm_dem_listy neighs[BMM_PART_MAX];
};

struct bmm_dem_link {
  // Link *to*.
  size_t i;
  // Spring rest position.
  double x0;
};

struct bmm_dem_listl {
  size_t n;
  struct bmm_dem_link linkl[BMM_GROUP_MAX];
};

struct bmm_dem_buf {
  size_t npart;
  struct bmm_dem_neigh neigh;
  struct bmm_dem_partc partcs[BMM_PART_MAX];
  struct bmm_dem_part parts[BMM_PART_MAX];
  // These are directed links *to* some particle.
  struct bmm_dem_listl links[BMM_PART_MAX];
  // TODO Do we want a asymmetric pair list or an symmetric list of lists?
  // Probably the latter, even though keeping it consistent takes work.
};

struct bmm_dem {
  struct bmm_dem_opts opts;
  gsl_rng* rng;
  enum bmm_dem_mode mode;
  size_t istep;
  double rext[2];
  // TODO Initial value system goes here.
  // TODO Force scheme goes here.
  void (* forcesch)(struct bmm_dem*);
  double (* fnormal)(struct bmm_dem const*, double const*, double const*);
  double (* ftang)(struct bmm_dem const*, double const*, double const*);
  // TODO Integration scheme goes here.
  void (* intsch)(struct bmm_dem*);
  // TODO Measurement system goes here.
  struct bmm_dem_est est;
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
  // Large object memory pool.
  union {
    // Which particles each cell owns, roughly.
    struct bmm_dem_list conts[BMM_CELL_MAX * BMM_CELL_MAX];
  } pool;
};

// TODO Combine list types.

// Index lists.

inline void bmm_dem_clear(struct bmm_dem_list* const list) {
  list->n = 0;
}

inline bool bmm_dem_push(struct bmm_dem_list* const list, size_t const x) {
  if (list->n >= sizeof list->i / sizeof *list->i)
    return false;

  list->i[list->n] = x;
  ++list->n;

  return true;
}

inline size_t bmm_dem_size(struct bmm_dem_list const* const list) {
  return list->n;
}

inline size_t bmm_dem_get(struct bmm_dem_list const* const list,
    size_t const i) {
  return list->i[i];
}

// TODO Some other kinds of lists.

inline void bmm_dem_cleary(struct bmm_dem_listy* const list) {
  list->n = 0;
}

inline bool bmm_dem_pushy(struct bmm_dem_listy* const list,
    size_t const x) {
  if (list->n >= sizeof list->thingy / sizeof *list->thingy)
    return false;

  list->thingy[list->n].i = x;
  list->thingy[list->n].x[0] = 0.0;
  list->thingy[list->n].x[1] = 0.0;
  ++list->n;

  return true;
}

inline size_t bmm_dem_sizey(struct bmm_dem_listy const* const list) {
  return list->n;
}

inline size_t bmm_dem_gety(struct bmm_dem_listy const* const list,
    size_t const i) {
  return list->thingy[i].i;
}

// TODO Yet another kinds of lists.

inline void bmm_dem_clearl(struct bmm_dem_listl* const list) {
  list->n = 0;
}

inline bool bmm_dem_pushl(struct bmm_dem_listl* const list,
    size_t const x) {
  if (list->n >= sizeof list->linkl / sizeof *list->linkl)
    return false;

  list->linkl[list->n].i = x;
  list->linkl[list->n].x0 = 0.0;
  ++list->n;

  return true;
}

inline size_t bmm_dem_sizel(struct bmm_dem_listl const* const list) {
  return list->n;
}

inline size_t bmm_dem_getl(struct bmm_dem_listl const* const list,
    size_t const i) {
  return list->linkl[i].i;
}

/// The call `bmm_dem_getbuf(dem)`
/// returns the active read and write buffer of the simulation `dem`
/// whether it is single-buffered or double-buffered.
__attribute__ ((__nonnull__))
inline struct bmm_dem_buf* bmm_dem_getbuf(struct bmm_dem* const dem) {
  return dem->dblbuf ? dem->data.bufs.active : &dem->data.buf;
}

/// The call `bmm_dem_getrbuf(dem)`
/// returns the active read buffer of the simulation `dem`
/// whether it is single-buffered or double-buffered.
__attribute__ ((__nonnull__))
inline struct bmm_dem_buf const* bmm_dem_getrbuf(
    struct bmm_dem const* const dem) {
  return dem->dblbuf ? dem->data.bufs.active : &dem->data.buf;
}

/// The call `bmm_dem_getwbuf(dem)`
/// returns the passive write buffer of the simulation `dem`
/// whether it is single-buffered or double-buffered.
__attribute__ ((__nonnull__))
inline struct bmm_dem_buf* bmm_dem_getwbuf(struct bmm_dem* const dem) {
  return dem->dblbuf ? dem->data.bufs.passive : &dem->data.buf;
}

/// The call `bmm_dem_swapbuf(dem)`
/// swaps the active and passive buffers of the simulation `dem`.
/// This is necessary for operations that are not structure-preserving.
__attribute__ ((__nonnull__))
inline void bmm_dem_swapbuf(struct bmm_dem* const dem) {
  if (dem->dblbuf) {
    struct bmm_dem_buf* const buf = dem->data.bufs.passive;
    dem->data.bufs.passive = dem->data.bufs.active;
    dem->data.bufs.active = buf;
  }
}

__attribute__ ((__nonnull__))
void bmm_dem_opts_def(struct bmm_dem_opts*);

__attribute__ ((__nonnull__))
void bmm_dem_part_def(struct bmm_dem_part*);

__attribute__ ((__nonnull__))
void bmm_dem_def(struct bmm_dem*, struct bmm_dem_opts const*);

__attribute__ ((__nonnull__))
double bmm_dem_ekinetic(struct bmm_dem const*);

__attribute__ ((__nonnull__))
double bmm_dem_pvector(struct bmm_dem const*);

__attribute__ ((__nonnull__))
double bmm_dem_pscalar(struct bmm_dem const*);

__attribute__ ((__nonnull__))
double bmm_dem_cor(struct bmm_dem const*);

__attribute__ ((__nonnull__))
bool bmm_dem_step(struct bmm_dem*);

__attribute__ ((__nonnull__))
bool bmm_dem_comm(struct bmm_dem*);

__attribute__ ((__nonnull__))
bool bmm_dem_run(struct bmm_dem*);

__attribute__ ((__nonnull__))
bool bmm_dem_run_with(struct bmm_dem_opts const*);

#endif
