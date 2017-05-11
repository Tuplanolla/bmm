#ifndef BMM_DEM_H
/// Discrete element method with some assumptions.
#define BMM_DEM_H

#include <gsl/gsl_rng.h>
#include <stdbool.h>
#include <stddef.h>

#include "conf.h"
#include "cpp.h"
#include "ext.h"
#include "io.h"
#include "msg.h"

enum bmm_dem_init {
  BMM_DEM_INIT_TRIAL,
  BMM_DEM_INIT_CUBIC,
  BMM_DEM_INIT_POISSOND
};

enum bmm_dem_integ {
  BMM_DEM_INTEG_EULER,
  BMM_DEM_INTEG_GEAR
};

enum bmm_dem_fnorm {
  BMM_DEM_FNORM_DASHPOT,
  BMM_DEM_FNORM_VISCOEL
};

enum bmm_dem_ftang {
  BMM_DEM_FTANG_HW,
  BMM_DEM_FTANG_CS
};

enum bmm_dem_role {
  BMM_DEM_ROLE_FREE,
  BMM_DEM_ROLE_FIXED,
  BMM_DEM_ROLE_DRIVEN
};

enum bmm_dem_mode {
  BMM_DEM_MODE_BEGIN,
  BMM_DEM_MODE_SEDIMENT,
  BMM_DEM_MODE_LINK,
  BMM_DEM_MODE_BREAK,
  BMM_DEM_MODE_ACCEL,
  BMM_DEM_MODE_CRUNCH
};

struct bmm_dem_opts {
  /// Initialization scheme.
  enum bmm_dem_init init;
  /// Integration scheme.
  enum bmm_dem_integ integ;
  /// Normal force scheme.
  enum bmm_dem_fnorm fnorm;
  /// Tangential force scheme.
  enum bmm_dem_ftang ftang;
  /// Bounding box.
  struct {
    /// Extents.
    double r[BMM_NDIM];
    /// Periodicities.
    bool per[BMM_NDIM];
  } box;
  /// Ambient properties.
  struct {
    /// Fictitious dynamic viscosity $\\eta$
    /// used for calculating the drag force $F = -b v$,
    /// where the drag coefficient $b = 3 \\twopi \\eta r$ and
    /// the Stokes radius $r$ is equal to particle radius.
    double eta;
  } ambient;
  /// Particles.
  struct {
    /// Young's modulus $Y$.
    double y;
    /// Dashpot elasticity $\\gamma$.
    double gamma;
    /// Particle sizes expressed as the support of the uniform distribution.
    double rnew[2];
  } part;
  /// Links between particles.
  struct {
    /// Link length creation factor.
    double ccrlink;
    /// Link length shrink factor.
    double cshlink;
    /// Tensile spring constant.
    double ktens;
    /// Shear spring constant.
    double kshear;
  } link;
  /// Script to follow.
  struct {
    /// Number of stages.
    size_t n;
    /// Modes and their timings for each stage.
    struct {
      /// Timespan.
      double tspan;
      /// Time step.
      double dt;
      /// Functionality.
      enum bmm_dem_mode mode;
      /// Parameters.
      union {
        /// For `BMM_DEM_MODE_CRUNCH`.
        struct {
          /// Driving velocity target.
          double v[BMM_NDIM];
          /// Force increment.
          double fadjust;
        } crunch;
        /// For `BMM_DEM_MODE_BREAK`.
        struct {
          /// Pull-apart vector.
          double xgap[BMM_NDIM];
        } smash;
        /// For `BMM_DEM_MODE_SEDIMENT`.
        struct {
          /// Cohesive force.
          double fcohes;
        } sediment;
      } params;
    } stage[BMM_NSTAGE];
  } script;
  /// Communications.
  struct {
    /// Time step.
    double dt;
    /// Send this.
    bool flip;
    /// Send that.
    bool flop;
    /// Send these.
    bool flap;
    /// Send those.
    bool flup;
  } comm;
  /// Neighbor cache tuning.
  struct {
    /// Number of cells for each dimension.
    size_t ncell[BMM_NDIM];
    /// Maximum distance for qualifying as a neighbor.
    double rcutoff;
  } cache;
};

struct bmm_dem {
  struct bmm_dem_opts opts;
  /// Random number generator state.
  gsl_rng* rng;
  /// Particle pair interactions.
  struct {
    /// Compressions.
    double xi[BMM_NPART][BMM_NPART];
    /// Compression velocities.
    double xip[BMM_NPART][BMM_NPART];
    /// Compression accelerations.
    double xipp[BMM_NPART][BMM_NPART];
  } pair;
  /// Particles.
  struct {
    /// Number of particles.
    size_t n;
    /// Next unused label.
    size_t lnew;
    /// Labels.
    size_t l[BMM_NPART];
    /// Roles.
    enum bmm_dem_role role[BMM_NPART];
    /// Radii.
    double r[BMM_NPART];
    /// Masses.
    double m[BMM_NPART];
    /// Moments of inertia.
    double j[BMM_NPART];
    /// Step.
    double i;
    /// Time.
    double t;
    /// Positions.
    double x[BMM_NPART][BMM_NDIM];
    /// Velocities.
    double v[BMM_NPART][BMM_NDIM];
    /// Accelerations.
    double a[BMM_NPART][BMM_NDIM];
    /// Angles.
    double phi[BMM_NPART];
    /// Angular velocities.
    double omega[BMM_NPART];
    /// Angular accelerations.
    double alpha[BMM_NPART];
    /// Forces.
    double f[BMM_NPART][BMM_NDIM];
    /// Torques.
    double tau[BMM_NPART];
  } part;
  /// Links between particles.
  struct {
    /// Number of links.
    size_t n;
    /// Index pairs, each in ascending order.
    size_t i[BMM_NLINK][2];
    /// Rest lengths for springs and beams.
    double rrest[BMM_NLINK];
    /// Rest angles for beams.
    double phirest[BMM_NLINK][2];
    /// Limit force for tensile stress induced breaking.
    double ftens[BMM_NLINK];
    /// Limit force for shear stress induced breaking.
    double fshear[BMM_NLINK];
  } link;
  /// Script state.
  struct {
    /// Current stage.
    size_t i;
    /// Previous transition time.
    double tprev;
    /// Next transition time.
    double tnext;
    /// State of the current mode.
    union {
      /// Total driving force for `BMM_DEM_MODE_CRUNCH`.
      double f[BMM_NDIM];
    } state;
  } script;
  /// Communications.
  struct {
    /// Current message.
    size_t i;
    /// Previous message time.
    double tprev;
    /// Next message time.
    double tnext;
  } comm;
  /// Neighbor cache.
  /// This is only used for performance optimization.
  struct {
    /// Time of previous update.
    double tprev;
    /// Time of next update.
    double tnext;
    /// Which particles were previously in each neighbor cell.
    struct {
      /// Number of particles.
      size_t n;
      /// Particle indices.
      size_t i[BMM_NGROUP];
    } part[BMM_POW(BMM_NCELL, BMM_NDIM)];
    /// Which neighbors each particle previously had.
    /// This only covers half of the Moore neighborhood of a particle.
    struct {
      /// Number of neighbors.
      size_t n;
      /// Neighbor indices.
      size_t i[BMM_NGROUP * (BMM_POW(3, BMM_NDIM) - 1)];
    } neigh[BMM_NPART];
  } cache;
};

// TODO Combine list types.

// Index lists.

/*
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
*/

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

// TODO These are questionable to expose.

size_t bmm_dem_sniff_size(struct bmm_dem const* const dem,
    enum bmm_msg_num const num);

bool bmm_dem_puts_stuff(struct bmm_dem const* const dem,
    enum bmm_msg_num const num);

bool bmm_dem_puts(struct bmm_dem const* const dem,
    enum bmm_msg_num const num);

#endif
