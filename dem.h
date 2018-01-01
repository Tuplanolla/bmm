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

/// Special particle properties.
enum bmm_dem_role {
  /// Free particle.
  BMM_DEM_ROLE_FREE,
  /// Fixed particle.
  BMM_DEM_ROLE_FIXED,
  /// Driven particle.
  BMM_DEM_ROLE_DRIVEN
};

/// Cache policies.
enum bmm_dem_cache {
  BMM_DEM_CACHE_NONE,
  /// Neighbor cell caching.
  BMM_DEM_CACHE_NEIGH
};

/// Contact types for pairs of particles.
enum bmm_dem_ct {
  /// Weak contact that may develop into no contact.
  BMM_DEM_CT_WEAK,
  /// Strong contact that may develop into weak or no contact.
  BMM_DEM_CT_STRONG,
  /// Number of contact types.
  BMM_NCT
};

/// Endpoints for pairs of particles.
enum bmm_dem_end {
  /// Back end.
  BMM_DEM_END_TAIL,
  /// Front end.
  BMM_DEM_END_HEAD,
  /// Number of ends.
  BMM_NEND
};

/// Integration schemes.
enum bmm_dem_integ {
  /// Forward Euler scheme.
  BMM_DEM_INTEG_EULER,
  /// Forward scheme using Taylor expansions.
  BMM_DEM_INTEG_TAYLOR,
  /// Velocity Verlet (not leapfrog) scheme.
  BMM_DEM_INTEG_VELVET,
  /// Beeman (modified) scheme.
  BMM_DEM_INTEG_BEEMAN,
  /// Kuraev (modified) scheme.
  BMM_DEM_INTEG_KURAEV
};

/// External force schemes.
enum bmm_dem_ext {
  BMM_DEM_EXT_NONE,
  /// Harmonic force field.
  BMM_DEM_EXT_HARM,
  /// Gravitational acceleration.
  BMM_DEM_EXT_GRAVY,
  /// Driving force.
  BMM_DEM_EXT_DRIVE
};

/// Ambient force schemes.
enum bmm_dem_amb {
  BMM_DEM_AMB_NONE,
  /// Free stream creeping flow.
  BMM_DEM_AMB_FAXEN,
  /// Laminar channel flow.
  BMM_DEM_AMB_LAMCHAN,
  /// Turbulent channel flow.
  BMM_DEM_AMB_TURBCHAN
};

/// Normal force schemes.
enum bmm_dem_norm {
  BMM_DEM_NORM_NONE,
  /// Ideal restitution model.
  BMM_DEM_NORM_IDEAL,
  /// Linear dashpot model named after Kelvin and Voigt.
  BMM_DEM_NORM_KV,
  /// Elastic model by Hertz.
  BMM_DEM_NORM_HERTZ,
  /// Viscoelastic model by Brilliantov, Spahn, Hertzsch and Poschel.
  BMM_DEM_NORM_BSHP,
  /// Piecewise model by Walton and Braun.
  BMM_DEM_NORM_WB
};

/// Tangential force schemes.
enum bmm_dem_tang {
  BMM_DEM_TANG_NONE,
  /// Linear model by Coulomb.
  // TODO This is `\\mu F` only.
  BMM_DEM_TANG_COULOMB,
  /// Spring model by Cundall and Strack.
  BMM_DEM_TANG_CS,
  /// Power law model by Haff and Werner.
  BMM_DEM_TANG_HW,
  /// Microscopic asperity model by Brilliantov, Spahn, Hertzsch and Poschel.
  // TODO This is strange.
  BMM_DEM_TANG_BSHP,
  /// Piecewise model by Walton and Braun.
  BMM_DEM_TANG_WB,
  /// Ideal beam model.
  BMM_DEM_TANG_BEAM,
  /// Bending beam model by Euler and Bernoulli.
  BMM_DEM_TANG_EB,
  /// Bending and shear beam model.
  BMM_DEM_TANG_SHEAR,
  /// Bending and rotary beam model by Rayleigh.
  BMM_DEM_TANG_RAYLEIGH,
  /// Bending, shear and rotary beam model by Timoshenko.
  BMM_DEM_TANG_TIMO
};

/// Bonding criteria.
enum bmm_dem_bond {
  BMM_DEM_BOND_NONE
};

/// Yield criteria.
enum bmm_dem_yield {
  BMM_DEM_YIELD_NONE,
  /// Maximum normal stress criterion by Rankine.
  BMM_DEM_YIELD_RANKINE,
  /// Maximum normal strain criterion by Venant.
  BMM_DEM_YIELD_VENANT,
  /// Maximum shear stress criterion by Tresca.
  BMM_DEM_YIELD_TRESCA,
  /// Maximum distortion energy criterion by von Mises.
  BMM_DEM_YIELD_VONMISES,
  /// Mohr--Coulomb criterion.
  BMM_DEM_YIELD_MC,
  /// Ellipse criterion by Zhang and Eckert.
  BMM_DEM_YIELD_ZE,
  /// Pressure criterion by Drucker and Prager.
  BMM_DEM_YIELD_DP,
  /// Extended pressure criterion by Bresler and Pister.
  BMM_DEM_YIELD_BP,
  /// Combined criterion by Willam and Warnke.
  BMM_DEM_YIELD_WW
};

enum bmm_dem_mode {
  /// Do nothing.
  BMM_DEM_MODE_IDLE,
  /// Adjust the particle density.
  BMM_DEM_MODE_SET_DENSITY,
  /// Create a fixed number of particles, sparse in the y-direction.
  BMM_DEM_MODE_CREATE_GAS,
  /// Test systems.
  BMM_DEM_MODE_CREATE_COUPLE,
  BMM_DEM_MODE_CREATE_TRIPLET,
  BMM_DEM_MODE_CREATE_HC,
  BMM_DEM_MODE_CREATE_HEX,
  BMM_DEM_MODE_CREATE_PILE,
  BMM_DEM_MODE_CREATE_PLANE,
  BMM_DEM_MODE_CREATE_BEAM,
  /// Draw particles towards a harmonic force field zero along the x-axis.
  BMM_DEM_MODE_SEDIMENT,
  /// Draw particles toward infinity along the y-axis.
  BMM_DEM_MODE_GRAVY,
  /// Remove particles outside the bounding box.
  BMM_DEM_MODE_CLIP,
  // TODO Choose contact strengths from the Weibull distribution.
  /// Link nearby particles together.
  BMM_DEM_MODE_LINK,
  /// Break the contacts that cross the zero along the x-axis.
  BMM_DEM_MODE_FAULT,
  /// Pull the material apart in the y-direction.
  BMM_DEM_MODE_SEPARATE,
  // TODO Also "fix" the top during this stage.
  /// Fix the bottom.
  BMM_DEM_MODE_GLUE,
  /// Start forcing the top in the x-direction.
  BMM_DEM_MODE_CRUNCH,
  /// Begin measurements.
  BMM_DEM_MODE_MEASURE,
  /// Debug utilities.
  BMM_DEM_MODE_PRESET0,
  BMM_DEM_MODE_PRESET1,
  BMM_DEM_MODE_PRESET2
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
  /// Fault.
  struct {
    /// Fault jags.
    size_t njag;
    /// Fault jag height.
    double hjag;
  } fault;
  /// Timekeeping.
  struct {
    /// Stabilization frequency (frame rule).
    size_t istab;
  } time;
  /// Particles.
  struct {
    /// True mass density.
    double rho;
    /// Tensile Young's modulus.
    double ytens;
    /// Compressive Young's modulus.
    double ycomp;
    /// Poisson's ratio.
    double nu;
    /// Particle sizes expressed as the width of the uniform distribution.
    double rnew[2];
    /// Link strengths expressed as the width of the uniform distribution.
    double strnew[2];
  } part;
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
      /// For `BMM_DEM_MODE_GLUE`.
      struct {
        /// Thickness of fixed layer.
        double nlayer;
      } glue;
      /// For `BMM_DEM_MODE_CRUNCH`.
      struct {
        /// Thickness of driven layer.
        double nlayer;
        /// Driving velocity target.
        double v;
        /// Force increment.
        double fadjust[BMM_NDIM];
        /// Pressure target.
        double p;
      } crunch;
      /// For `BMM_DEM_MODE_FAULT`.
      struct {
        /// Fault shape (via indication).
        bool (*ind)(double const *, void *);
      } fault;
      /// For `BMM_DEM_MODE_SEPARATE`.
      struct {
        /// Pull-apart vector.
        double xgap[BMM_NDIM];
      } separate;
      /// For `BMM_DEM_MODE_SEDIMENT`.
      struct {
        /// Strength of cohesive force.
        double kcohes;
      } sediment;
      /// For `BMM_DEM_MODE_GRAVY`.
      struct {
        /// Gravitational acceleration.
        double g;
      } gravy;
      /// For `BMM_DEM_MODE_PRESET0` and `BMM_DEM_MODE_PRESET1`.
      struct {
        // Whether to use HW or CS.
        bool statfric;
        // As defined.
        double eta;
        double eta2;
        double a;
        double mu;
        double kt;
        double gammat;
        double kn;
        double gamman;
        double barkt;
        double bargammat;
        double sigmacrit;
        double taucrit;
        double sigmacritt;
        double taucritt;
      } preset;
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
    double dcutoff;
  } cache;
};

struct bmm_dem_pair {
  /// Contacts between particles.
  struct {
    /// Source mappings.
    ///
    /// Theoretically contacts are symmetric and reflexive,
    /// but here they are represented asymmetrically and nonreflexively
    /// such that the source index is always smaller than the target index.
    /// Theoretically contacts may also have unbounded effective distances,
    /// but here they may be limited by practical constraints
    /// such as `opts.cache.dcutoff`.
    struct {
      /// Number of targets.
      size_t n;
      /// Target indices.
      size_t itgt[BMM_MCONTACT];
      /// Rest distances.
      double drest[BMM_MCONTACT];
      /// Rest angle.
      double psirest[BMM_MLINK][BMM_NEND];
      /// Strength scale factor.
      double strength[BMM_MCONTACT];
      /// Fatigue time (very small).
      size_t tfat[BMM_MCONTACT];
    } src[BMM_MPART];
  } cont;
  /// Only contact or more.
  bool cohesive;
  /// Normal forces.
  struct {
    /// Force scheme.
    enum bmm_dem_norm tag;
    /// Parameters.
    union {
      /// For `BMM_DEM_NORM_KV`.
      struct {
        /// Dashpot elasticity.
        double k;
        /// Dashpot viscosity.
        double gamma;
      } dashpot;
      /// For `BMM_DEM_NORM_BSHP`.
      struct {
        /// Dissipative constant.
        double a;
      } viscoel;
    } params;
  } norm;
  /// Tangential forces.
  struct {
    /// Force scheme.
    enum bmm_dem_tang tag;
    /// Parameters.
    union {
      /// For `BMM_DEM_TANG_HW`.
      struct {
        /// Haff--Werner elasticity.
        double gamma;
        /// Coulomb friction parameter.
        double mu;
      } hw;
      /// For `BMM_DEM_TANG_CS`.
      struct {
        /// Cundall--Strack elasticity.
        double k;
        /// Coulomb friction parameter.
        double mu;
      } cs;
      /// For `BMM_DEM_TANG_BEAM`.
      struct {
        /// Stiffness.
        double k;
        /// Damping.
        double dk;
      } beam;
    } params;
  } tang;
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
        /// Old accelerations.
        double ao[BMM_MPART][BMM_NDIM];
        /// Old angular accelerations.
        double alphao[BMM_MPART];
      } velvet;
      /// For `BMM_DEM_INTEG_BEEMAN`.
      struct {
        /// Old displacements.
        double xo[BMM_MPART][BMM_NDIM];
        /// Old velocities.
        double vo[BMM_MPART][BMM_NDIM];
        /// Old accelerations.
        double ao[BMM_MPART][BMM_NDIM];
        /// Very old accelerations.
        double aoo[BMM_MPART][BMM_NDIM];
        /// Old angular displacements.
        double phio[BMM_MPART];
        /// Old angular velocities.
        double omegao[BMM_MPART];
        /// Old angular accelerations.
        double alphao[BMM_MPART];
        /// Very old angular accelerations.
        double alphaoo[BMM_MPART];
      } beeman;
    } params;
  } integ;
  /// External forces.
  struct {
    /// Force scheme.
    enum bmm_dem_ext tag;
    /// Parameters.
    union {
      /// For `BMM_DEM_EXT_HARM`.
      struct {
        /// Harmonic constant.
        double k;
        // TODO Add directional parameters here.
      } harm;
      /// For `BMM_DEM_EXT_GRAVY`.
      struct {
        /// Acceleration.
        double g;
      } gravy;
    } params;
  } ext;
  /// Ambient forces.
  struct {
    /// Force scheme.
    enum bmm_dem_amb tag;
    /// Parameters.
    union {
      /// For `BMM_DEM_AMB_CREEPING`.
      struct {
        /// Dynamic viscosity of compressible solution.
        double eta;
      } creeping;
    } params;
  } amb;
  /// Contacts.
  struct bmm_dem_pair pair[BMM_NCT];
  /// Weak-to-strong bonding criteria.
  struct {
    /// Bonding criterion.
    enum bmm_dem_bond tag;
    /// Link length creation factor.
    double ccrcont;
  } bond;
  /// Strong-to-weak yield criteria.
  struct {
    /// Yield criterion.
    enum bmm_dem_yield tag;
    /// Parameters.
    union {
      /// For `BMM_DEM_YIELD_RANKINE`.
      struct {
        /// Maximum normal stress.
        double sigmacrit;
      } rankine;
      /// For `BMM_DEM_YIELD_TRESCA`.
      struct {
        /// Maximum shear stress.
        double taucrit;
      } tresca;
      /// For `BMM_DEM_YIELD_ZE`.
      struct {
        /// Maximum normal stress.
        double sigmacrit;
        /// Maximum shear stress.
        double taucrit;
        /// Maximum normal stress (tensile).
        double sigmacritt;
        /// Maximum shear stress (tensile).
        double taucritt;
      } ze;
    } params;
  } yield;
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
        /// Number of fixed particles.
        size_t nfix;
        /// Number of driven particles.
        size_t ndrive;
        /// Total driving force.
        double fdrive[BMM_NDIM];
      } crunch;
    } state;
  } script;
  /// Communications.
  struct {
    /// Previous message time.
    double tprev;
  } comm;
  /// Estimator cache.
  /// This is only used for programmer laziness.
  struct {
    /// Ambient energy loss.
    double eambdis;
    /// The `_d` prefix implies (path-independent) direct measurements;
    /// they are more accurate than integrated measurements,
    /// with the exception that bonding and yielding makes them bogus.
    /// External potential energy.
    double epotext_d;
    /// Linear kinetic energy.
    double eklin_d;
    /// Rotational kinetic energy.
    double ekrot_d;
    /// Weak contact energy.
    double ewcont_d;
    /// Strong contact energy.
    double escont_d;
    /// Weak contact energy.
    double ewcont;
    /// Strong contact energy.
    double escont;
    /// Work done by cohesive pressure.
    double edrivnorm;
    /// Work done by driving force.
    double edrivtang;
    /// Energy forced into bonding.
    double ebond;
    /// Energy dissipated in yielding.
    double eyieldis;
    /// Energy dissipated in weak contact.
    double ewcontdis;
    /// Energy dissipated in strong contact.
    double escontdis;
    /// Resisting force.
    double fback[BMM_NDIM];
    /// Fixing velocity.
    double vfix[BMM_NDIM];
    /// Driving velocity.
    double vdriv[BMM_NDIM];
    /// Observed packing fraction.
    double chi;
    /// Effective macroscopic friction factor.
    double mueff;
    /// Effective macroscopic friction factor (from feedback).
    double mueffb;
  } est;
  /// Neighbor cache.
  /// This is only used for performance optimization.
  struct {
    /// Caching scheme.
    enum bmm_dem_cache tag;
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

__attribute__ ((__nonnull__))
void bmm_dem_opts_set_rnew(struct bmm_dem_opts *, double const *);

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
