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

// TODO Not sure if good idea.
enum bmm_dem_step {
  BMM_DEM_STEP_ERROR,
  BMM_DEM_STEP_STOP,
  BMM_DEM_STEP_CONT
};

enum bmm_dem_init {
  BMM_DEM_INIT_NONE,
  BMM_DEM_INIT_TRIAL,
  BMM_DEM_INIT_CUBIC,
  BMM_DEM_INIT_POISSOND
};

enum bmm_dem_integ {
  BMM_DEM_INTEG_EULER,
  BMM_DEM_INTEG_GEAR
};

enum bmm_dem_caching {
  BMM_DEM_CACHING_NONE,
  BMM_DEM_CACHING_NEIGH
};

enum bmm_dem_famb {
  BMM_DEM_FAMB_NONE,
  BMM_DEM_FAMB_CREEPING
};

enum bmm_dem_fnorm {
  BMM_DEM_FNORM_NONE,
  BMM_DEM_FNORM_DASHPOT,
  BMM_DEM_FNORM_VISCOEL
};

enum bmm_dem_ftang {
  BMM_DEM_FTANG_NONE,
  BMM_DEM_FTANG_HW,
  BMM_DEM_FTANG_CS
};

enum bmm_dem_flink {
  BMM_DEM_FLINK_NONE,
  BMM_DEM_FLINK_SPRING,
  BMM_DEM_FLINK_BEAM
};

enum bmm_dem_role {
  BMM_DEM_ROLE_FREE,
  BMM_DEM_ROLE_FIXED,
  BMM_DEM_ROLE_DRIVEN
};

enum bmm_dem_mode {
  BMM_DEM_MODE_IDLE,
  BMM_DEM_MODE_BEGIN,
  BMM_DEM_MODE_SEDIMENT,
  BMM_DEM_MODE_LINK,
  BMM_DEM_MODE_SMASH,
  BMM_DEM_MODE_ACCEL,
  BMM_DEM_MODE_CRUNCH
};

struct bmm_dem_opts {
  /// Report statistics.
  bool verbose;
  /// Initialization scheme.
  enum bmm_dem_init init;
  /// Integration scheme.
  enum bmm_dem_integ integ;
  /// Caching scheme.
  enum bmm_dem_caching caching;
  /// Ambient force scheme.
  enum bmm_dem_famb famb;
  /// Normal force scheme.
  enum bmm_dem_fnorm fnorm;
  /// Tangential force scheme.
  enum bmm_dem_ftang ftang;
  /// Link force scheme.
  enum bmm_dem_flink flink;
  /// Bounding box.
  struct {
    /// Extents.
    double x[BMM_NDIM];
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
  /// Normal forces.
  struct {
    /// Parameters.
    union {
      /// For `BMM_DEM_FNORM_DASHPOT`.
      struct {
        /// Dashpot elasticity.
        double gamma;
      } dashpot;
    } params;
  } norm;
  /// Tangential forces.
  struct {
    /// Parameters.
    union {
      /// For `BMM_DEM_FTANG_HW`.
      struct {
        /// Haff--Werner elasticity.
        double gamma;
        /// Coulomb friction parameter.
        double mu;
      } hw;
    } params;
  } tang;
  /// Particles.
  struct {
    /// Young's modulus.
    double y;
    /// Particle sizes expressed as the support of the uniform distribution.
    double rnew[2];
  } part;
  /// Links between particles.
  struct {
    /// Link length creation factor.
    double ccrlink;
    /// Link length expansion factor.
    double cshlink;
    /// Tensile spring constant.
    double ktens;
    /// Shear spring constant.
    double kshear;
    /// Limit length factors for tensile stress induced breaking
    /// expressed as the support of the uniform distribution.
    double crlim[2];
    /// Limit angle factors for shear stress induced breaking
    /// expressed as the support of the uniform distribution.
    double cphilim[2];
  } link;
  /// Script to follow.
  struct {
    /// Number of stages.
    size_t n;
    /// Timespans.
    double tspan[BMM_MSTAGE];
    /// Time steps.
    double dt[BMM_MSTAGE];
    /// Functionalities.
    enum bmm_dem_mode mode[BMM_MSTAGE];
    /// Parameters.
    union {
      /// For `BMM_DEM_MODE_CRUNCH`.
      struct {
        /// Driving velocity target.
        double v[BMM_NDIM];
        /// Force increment.
        double fadjust;
      } crunch;
      /// For `BMM_DEM_MODE_SMASH`.
      struct {
        /// Pull-apart vector.
        double xgap[BMM_NDIM];
      } smash;
      /// For `BMM_DEM_MODE_SEDIMENT`.
      struct {
        /// Cohesive force.
        double fcohes;
      } sediment;
    } params[BMM_MSTAGE];
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
    /// Number of neighbor cells for each dimension.
    /// There are always at least $3^d$ neighbor cells,
    /// because those outside the bounding extend to infinity.
    size_t ncell[BMM_NDIM];
    /// Maximum distance for qualifying as a neighbor.
    double rcutoff;
  } cache;
};

struct bmm_dem {
  struct bmm_dem_opts opts;
  /// Random number generator state.
  gsl_rng* rng;
  /// Timekeeping.
  struct {
    /// Step.
    size_t istep;
    /// Stabilization frequency (frame rule).
    size_t istab;
    /// Time.
    double t;
  } time;
  /// Particles.
  struct {
    /// Number of particles.
    size_t n;
    /// Next unused label.
    size_t lnew;
    /// Labels.
    size_t l[BMM_MPART];
    /// Roles.
    enum bmm_dem_role role[BMM_MPART];
    /// Radii.
    double r[BMM_MPART];
    /// Masses.
    double m[BMM_MPART];
    /// Moments of inertia.
    double j[BMM_MPART];
    /// Positions.
    double x[BMM_MPART][BMM_NDIM];
    /// Velocities.
    double v[BMM_MPART][BMM_NDIM];
    /// Accelerations.
    double a[BMM_MPART][BMM_NDIM];
    /// Angles.
    double phi[BMM_MPART];
    /// Angular velocities.
    double omega[BMM_MPART];
    /// Angular accelerations.
    double alpha[BMM_MPART];
    /// Forces.
    double f[BMM_MPART][BMM_NDIM];
    /// Torques.
    double tau[BMM_MPART];
  } part;
  /// Links between particles.
  struct {
    /// Which other particles each particle is linked to.
    struct {
      /// Number of links from this particle.
      size_t n;
      /// Target particle indices.
      size_t i[BMM_MLINK];
      /// Rest lengths for springs and beams.
      double rrest[BMM_MLINK];
      /// Rest angles for beams.
      double phirest[BMM_MLINK][2];
      /// Limit length for tensile stress induced breaking.
      double rlim[BMM_MLINK];
      /// Limit angle for shear stress induced breaking.
      double philim[BMM_MLINK];
    } part[BMM_MPART];
  } link;
  /// Script state.
  struct {
    /// Current stage (may be one past the end to signal the end).
    size_t i;
    /// Previous transition time.
    double tprev;
    /// Transition times away from states.
    double ttrans[BMM_MSTAGE];
    /// Transition time offsets (positive means transition was late).
    double toff[BMM_MSTAGE];
    /// State of the current mode.
    union {
      /// For `BMM_DEM_MODE_CRUNCH`.
      struct {
        /// Total driving force.
        double f[BMM_NDIM];
      } crunch;
    } state;
  } script;
  /// Communications.
  struct {
    /// Next unused message.
    size_t lnew;
    /// Previous message time.
    double tprev;
    /// Next message time.
    double tnext;
  } comm;
  /// Neighbor cache.
  /// This is only used for performance optimization.
  struct {
    /// Current revision.
    size_t i;
    /// Time of previous partial update.
    double tpart;
    /// Time of previous full update.
    double tprev;
    /// Time of next full update.
    double tnext;
    /// Which neighbor cell each particle was in previously.
    size_t cell[BMM_MPART][BMM_NDIM];
    /// Which particles were previously in each neighbor cell.
    struct {
      /// Number of particles.
      size_t n;
      /// Particle indices.
      size_t i[BMM_MGROUP];
    } part[BMM_POW(BMM_MCELL, BMM_NDIM)];
    /// Which neighbors each particle previously had.
    /// This only covers half of the Moore neighborhood of a particle.
    struct {
      /// Number of neighbors.
      size_t n;
      /// Neighbor indices.
      size_t i[BMM_MGROUP * (BMM_POW(3, BMM_NDIM) / 2)];
    } neigh[BMM_MPART];
  } cache;
};

/// The call `bmm_dem_ijcell(pijcell, dem, ipart)`
/// writes the neighbor cell indices of the particle with the index `ipart`
/// into the index vector `pijcell`.
__attribute__ ((__nonnull__))
void bmm_dem_ijcell(size_t*, struct bmm_dem const*, size_t);

/// The call `bmm_dem_inspart(dem, r, m)`
/// places a new particle with radius `r` and mass `m`
/// in the origin at rest and
/// returns the index of the new particle.
/// Otherwise `BMM_MPART` is returned.
__attribute__ ((__nonnull__))
size_t bmm_dem_inspart(struct bmm_dem*, double, double);

/// The call `bmm_dem_delpart(dem, ipart)`
/// removes the particle with the index `ipart`.
/// Note that the index may be immediately assigned to another particle,
/// so all index caches should be purged.
/// This operation may be slow due to index reassignment.
__attribute__ ((__nonnull__))
bool bmm_dem_delpart(struct bmm_dem*, size_t);

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
enum bmm_dem_step bmm_dem_step(struct bmm_dem*);

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
