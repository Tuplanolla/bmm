/// Discrete element method with some assumptions.

#ifndef BMM_DEM_H
#define BMM_DEM_H

#include <gsl/gsl_rng.h>
#include <stdbool.h>
#include <stddef.h>

#include "conf.h"
#include "cpp.h"
#include "geom2d.h"
#include "ext.h"
#include "fp.h"
#include "io.h"
#include "msg.h"

enum bmm_dem_integ {
  /// Forward Euler.
  BMM_DEM_INTEG_EULER,
  /// Taylor series.
  BMM_DEM_INTEG_TAYLOR,
  /// Velocity Verlet.
  BMM_DEM_INTEG_VELVET
};

enum bmm_dem_caching {
  BMM_DEM_CACHING_NONE,
  BMM_DEM_CACHING_NEIGH
};

enum bmm_dem_fext {
  BMM_DEM_FEXT_NONE,
  BMM_DEM_FEXT_ABS,
  BMM_DEM_FEXT_HARM
};

enum bmm_dem_famb {
  BMM_DEM_FAMB_NONE,
  BMM_DEM_FAMB_CREEPING,
  BMM_DEM_FAMB_QUAD,
  BMM_DEM_FAMB_CORR
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

enum bmm_dem_tau {
  BMM_DEM_TAU_SOFT,
  BMM_DEM_TAU_HARD,
  BMM_DEM_TAU_AVERAGE,
  BMM_DEM_TAU_HALFWAY
};

enum bmm_dem_fract {
  BMM_DEM_FRACT_NONE,
  BMM_DEM_FRACT_ELLIPSE
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
  /// Do nothing.
  BMM_DEM_MODE_IDLE,
  /// Create a fixed number of particles, sparse in the y-direction.
  BMM_DEM_MODE_CREATE_GAS,
  BMM_DEM_MODE_CREATE_COUPLE,
  BMM_DEM_MODE_CREATE_TRIPLET,
  BMM_DEM_MODE_CREATE_HC,
  BMM_DEM_MODE_CREATE_HEX,
  BMM_DEM_MODE_CREATE_PILE,
  BMM_DEM_MODE_CREATE_BEAM,
  /// Draw particles towards a harmonic force field zero along the x-axis.
  BMM_DEM_MODE_SEDIMENT,
  BMM_DEM_MODE_GRAVY,
  /// Remove particles outside the bounding box.
  BMM_DEM_MODE_CLIP,
  /// Link nearby particles together.
  BMM_DEM_MODE_LINK,
  // TODO Consider another way to do this:
  // use a partitioning indicator function when linking.
  /// Break the links that cross the zero along the x-axis.
  BMM_DEM_MODE_FAULT,
  /// Pull the material apart in the y-direction.
  BMM_DEM_MODE_SEPARATE,
  /// Fix the bottom and start forcing the top in the x-direction.
  BMM_DEM_MODE_CRUNCH,
  /// Begin measurements.
  BMM_DEM_MODE_MEASURE,
  /// Debug utilities.
  BMM_DEM_MODE_PRESET0,
  BMM_DEM_MODE_PRESET1
};

struct bmm_dem_opts {
  /// Report statistics.
  bool verbose;
  /// Floating-point exceptions.
  struct {
    /// Trap exceptions.
    bool enabled;
    /// Exception mask.
    int mask;
  } trap;
  /// Bounding box.
  struct {
    /// Extents.
    double x[BMM_NDIM];
    /// Periodicities.
    bool per[BMM_NDIM];
  } box;
  /// Timekeeping.
  struct {
    /// Stabilization frequency (frame rule).
    size_t istab;
  } time;
  /// Particles.
  struct {
    /// Mass density.
    double rho;
    /// Young's modulus.
    double y;
    /// Poisson's ratio.
    double nu;
    /// Particle sizes expressed as the width of the uniform distribution.
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
    /// Tensile spring constant derivative.
    double dktens;
    /// Shear spring constant.
    double kshear;
    /// Shear spring constant derivative.
    double dkshear;
    /// Limit length factors for tensile stress induced breaking
    /// expressed as the width of the uniform distribution.
    double crlim[2];
    /// Limit angle factors for shear stress induced breaking
    /// expressed as the width of the uniform distribution.
    double cphilim[2];
  } link;
  /// Script to follow.
  struct {
    /// Number of stages.
    size_t n;
    /// Functionalities.
    enum bmm_dem_mode mode[BMM_MSTAGE];
    /// Timespans.
    double tspan[BMM_MSTAGE];
    /// Time steps.
    double dt[BMM_MSTAGE];
    /// Parameters.
    union {
      /// For `BMM_DEM_MODE_CREATE`.
      struct {
        /// Target packing fraction.
        double eta;
        /// Spacing factor.
        double cspace;
      } create;
      /// For `BMM_DEM_MODE_CREATE_*`.
      struct {
        /// Target packing fraction.
        double eta;
        /// Thickness.
        double layers;
        /// Lengthness.
        double slices;
      } test;
      /// For `BMM_DEM_MODE_CRUNCH`.
      struct {
        /// Driving velocity target.
        double v[BMM_NDIM];
        /// Force increment.
        double fadjust;
      } crunch;
      /// For `BMM_DEM_MODE_FAULT`.
      struct {
        /// Pull-apart vector.
        double xgap[BMM_NDIM];
      } fault;
      /// For `BMM_DEM_MODE_SEDIMENT`.
      struct {
        /// Strength of cohesive force.
        double kcohes;
      } sediment;
      /// For `BMM_DEM_MODE_GRAVY`.
      struct {
        /// Gravitational force.
        double f;
      } gravy;
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
  gsl_rng *rng;
  /// Floating-point exceptions.
  struct {
    /// Original mask to restore.
    int remask;
  } trap;
  /// Predictor data.
  struct {
    /// Integration scheme.
    enum bmm_dem_integ tag;
    /// Parameters.
    union {
      /// For `BMM_DEM_INTEG_VELVET`.
      struct {
        /// Accelerations.
        double a[BMM_MPART][BMM_NDIM];
        /// Angular accelerations.
        double alpha[BMM_MPART];
      } velvet;
    } params;
  } integ;
  /// External forces.
  struct {
    /// External force scheme.
    enum bmm_dem_fext tag;
    /// Parameters.
    union {
      /// For `BMM_DEM_FEXT_ABS`.
      struct {
        /// Force.
        double fcohes;
      } abs;
      /// For `BMM_DEM_FEXT_HARM`.
      struct {
        /// Harmonic constant.
        double kcohes;
      } harm;
    } params;
  } ext;
  /// Ambient forces.
  struct {
    /// Ambient force scheme.
    enum bmm_dem_famb tag;
    /// Parameters.
    union {
      /// For `BMM_DEM_FAMB_CREEPING`.
      struct {
        /// Dynamic viscosity of compressible solution.
        double mu;
      } creeping;
    } params;
  } amb;
  /// Normal forces.
  struct {
    /// Normal force scheme.
    enum bmm_dem_fnorm tag;
    /// Parameters.
    union {
      /// For `BMM_DEM_FNORM_DASHPOT`.
      struct {
        /// Dashpot elasticity.
        double gamma;
      } dashpot;
      /// For `BMM_DEM_FNORM_VISCOEL`.
      struct {
        /// Dissipative constant.
        double a;
      } viscoel;
    } params;
  } norm;
  /// Tangential forces.
  struct {
    /// Tangential force scheme.
    enum bmm_dem_ftang tag;
    /// Parameters.
    union {
      /// For `BMM_DEM_FTANG_HW`.
      struct {
        /// Haff--Werner elasticity.
        double gamma;
        /// Coulomb friction parameter.
        double mu;
      } hw;
      /// For `BMM_DEM_FTANG_CS`.
      struct {
        /// Cundall--Strack elasticity.
        double kappa;
        /// Coulomb friction parameter.
        double mu;
        // TODO Pick a better data structure.
        /// Elongations.
        double zeta[BMM_MPART][BMM_MPART];
      } cs;
    } params;
  } tang;
  /// Torques.
  struct {
    /// Torque scheme.
    enum bmm_dem_tau tag;
  } tau;
  /// Timekeeping.
  struct {
    /// Time.
    double t;
    /// Step.
    size_t istep;
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
    /// Reduced moments of inertia.
    double jred[BMM_MPART];
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
  /// Link fracture.
  struct {
    /// Link force scheme.
    enum bmm_dem_fract tag;
  } fract;
  /// Links between particles.
  struct {
    /// Link force scheme.
    enum bmm_dem_flink tag;
    // TODO These four are unused for now.
    /// Link elasticity.
    double k;
    /// Link damping.
    double dk;
    /// Link angular elasticity.
    double kappa;
    /// Link angular damping.
    double dkappa;
    /// Which other particles each particle is linked to.
    struct {
      /// Number of links from this particle.
      size_t n;
      /// Target particle indices.
      size_t i[BMM_MLINK];
      /// Rest lengths for springs and beams.
      double rrest[BMM_MLINK];
      /// Rest angles for beams.
      double chirest[BMM_MLINK][2];
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
    /// Transition times (away from states).
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
    /// Previous message time.
    double tprev;
  } comm;
  /// Neighbor cache.
  /// This is only used for performance optimization.
  struct {
    /// Caching scheme.
    enum bmm_dem_caching tag;
    /// Freshness right now.
    bool stale;
    /// Current revision.
    size_t i;
    /// Time of previous partial update.
    double tpart;
    /// Time of previous full update.
    double tprev;
    /// Moments of inertia.
    double j[BMM_MPART];
    /// Previous positions.
    double x[BMM_MPART][BMM_NDIM];
    /// Which neighbor cell each particle was in previously.
    size_t ijcell[BMM_MPART][BMM_NDIM];
    /// Which neighbor cell each particle was in previously, with vengeance.
    size_t icell[BMM_MPART];
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
      size_t i[BMM_MGROUP * (BMM_POW(3, BMM_NDIM) / 2 + 1)];
    } neigh[BMM_MPART];
  } cache;
};

/// The call `bmm_dem_script_addstage(opts)`
/// tries to add an idle stage with zero duration to the script
/// in the simulation options `opts`.
/// If there is enough capacity and the operation is successful,
/// the index of the new stage is returned.
/// Otherwise `SIZE_MAX` is returned.
/// Note that zero duration means that the stage
/// takes exactly one frame during which no time elapses.
__attribute__ ((__nonnull__))
size_t bmm_dem_script_addstage(struct bmm_dem_opts *);

/// The call `bmm_dem_script_ongoing(dem)`
/// checks whether the simulation `dem` has not ended.
__attribute__ ((__nonnull__, __pure__))
bool bmm_dem_script_ongoing(struct bmm_dem const *);

/// The call `bmm_dem_script_trans(dem)`
/// transitions the simulation `dem` to the next state and
/// checks whether the simulation did not end while doing so.
/// Make sure the simulation has not ended prior to the call
/// by calling `bmm_dem_script_ongoing`.
__attribute__ ((__nonnull__))
bool bmm_dem_script_trans(struct bmm_dem *);

/// The call `bmm_dem_cache_build(dem)`
/// tries to cache everything that supports caching
/// in the simulation `dem`.
/// If there is enough capacity and the operation is successful,
/// `true` is returned.
/// Otherwise `false` is returned and
/// the cache is left in an undefined state.
__attribute__ ((__nonnull__))
bool bmm_dem_cache_build(struct bmm_dem *);

/// The call `bmm_dem_addpart(dem)`
/// tries to place a new particle with unit radius and unit mass
/// at rest at the origin
/// in the simulation `dem`.
/// If there is enough capacity and the operation is successful,
/// the index of the new particle is returned.
/// Otherwise `SIZE_MAX` is returned.
/// Note that placing new particles triggers rebuilding all caches
/// by calling `bmm_dem_cache_build`.
__attribute__ ((__nonnull__))
size_t bmm_dem_addpart(struct bmm_dem *);

/// The call `bmm_dem_rempart(dem, ipart)`
/// removes the particle with the index `ipart`
/// in the simulation `dem`.
/// Indices may be reassigned in the process.
/// Note that removing particles triggers rebuilding all caches
/// by calling `bmm_dem_cache_build`.
__attribute__ ((__nonnull__))
void bmm_dem_rempart(struct bmm_dem *, size_t);

/// The call `bmm_dem_force_pair(dem, ipart, jpart)`
/// calculates the forces between the particles `ipart` and `jpart`
/// in the simulation `dem`.
__attribute__ ((__nonnull__))
void bmm_dem_force_pair(struct bmm_dem *, size_t, size_t);

/// The call `bmm_dem_opts_def(opts)`
/// writes the default simulation options into `opts`.
__attribute__ ((__nonnull__))
void bmm_dem_opts_def(struct bmm_dem_opts *);

/// The call `bmm_dem_def(dem, opts)`
/// writes the default simulation state into `dem`
/// with the simulation options `opts`.
__attribute__ ((__nonnull__))
void bmm_dem_def(struct bmm_dem *, struct bmm_dem_opts const *);

double bmm_dem_est_ekin(struct bmm_dem const *);
double bmm_dem_est_cor(struct bmm_dem const *);

__attribute__ ((__nonnull__))
bool bmm_dem_step(struct bmm_dem *);

__attribute__ ((__nonnull__))
bool bmm_dem_comm(struct bmm_dem *);

/// The call `bmm_dem_stab(dem)`
/// stabilizes the simulation `dem`
/// by wrapping periodic values, renormalizing normal vectors and so on.
/// This should only affect the behavior of the simulation over long timespans.
__attribute__ ((__nonnull__))
void bmm_dem_stab(struct bmm_dem *);

/// The call `bmm_dem_report(dem)`
/// prints informal diagnostics for the simulation state `dem`.
__attribute__ ((__nonnull__))
bool bmm_dem_report(struct bmm_dem const *);

/// The call `bmm_dem_trap_on(dem)`
/// enables floating-point exception trapping
/// in the simulation `dem`.
__attribute__ ((__nonnull__))
bool bmm_dem_trap_on(struct bmm_dem *);

/// The call `bmm_dem_trap_off(dem)`
/// disables floating-point exception trapping and
/// restores the original state of the processor
/// in the simulation `dem`.
__attribute__ ((__nonnull__))
bool bmm_dem_trap_off(struct bmm_dem *);

/// The call `bmm_dem_run(dem)`
/// processes all incoming messages and
/// handles signals with the simulation state `dem`.
__attribute__ ((__nonnull__))
bool bmm_dem_run(struct bmm_dem *);

/// The call `bmm_dem_run_with(opts)`
/// processes all incoming messages and
/// handles signals with the simulation options `opts`.
__attribute__ ((__nonnull__))
bool bmm_dem_run_with(struct bmm_dem_opts const *);

// TODO These are questionable to expose.

size_t bmm_dem_sniff_size(struct bmm_dem const *const dem,
    enum bmm_msg_num const num);

bool bmm_dem_puts_stuff(struct bmm_dem const *const dem,
    enum bmm_msg_num const num);

bool bmm_dem_puts(struct bmm_dem const *const dem,
    enum bmm_msg_num const num);

/// The call `bmm_dem_cache_expired(dem)`
/// checks whether particles have moved enough to require rebuilding caches
/// in the simulation `dem`.
__attribute__ ((__nonnull__, __pure__))
bool bmm_dem_cache_expired(struct bmm_dem const *);

#endif
