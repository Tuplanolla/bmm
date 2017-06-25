#ifndef BMM_DEM_H
/// Discrete element method with some assumptions.
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
  BMM_DEM_MODE_CREATE,
  /// Draw particles towards a harmonic force field zero along the x-axis.
  BMM_DEM_MODE_SEDIMENT,
  /// Link nearby particles together.
  BMM_DEM_MODE_LINK,
  // TODO Consider another way to do this:
  // use an partitioning indicator function when linking.
  /// Break the links that cross the zero along the x-axis.
  BMM_DEM_MODE_FAULT,
  /// Pull the material apart in the y-direction.
  BMM_DEM_MODE_SEPARATE,
  /// Fix the bottom and start forcing the top in the x-direction.
  BMM_DEM_MODE_CRUNCH,
  /// Begin measurements.
  BMM_DEM_MODE_MEASURE
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
    /// Parameters.
    union {
      /// For `BMM_DEM_FAMB_CREEPING`.
      struct {
        /// Dynamic viscosity of compressible solution.
        double mu;
      } creeping;
    } params;
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
  /// Timekeeping.
  struct {
    /// Stabilization frequency (frame rule).
    size_t istab;
  } time;
  /// Particles.
  struct {
    /// Young's modulus.
    double y;
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
    /// Shear spring constant.
    double kshear;
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
    /// Timespans.
    double tspan[BMM_MSTAGE];
    /// Time steps.
    double dt[BMM_MSTAGE];
    /// Functionalities.
    enum bmm_dem_mode mode[BMM_MSTAGE];
    /// Parameters.
    union {
      /// For `BMM_DEM_MODE_CREATE`.
      struct {
        /// Target packing fraction.
        double eta;
      } create;
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
    size_t cell[BMM_MPART][BMM_NDIM];
    /// Which neighbor cell each particle was in previously, with vengeance.
    size_t hurr[BMM_MPART];
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

/// The call `bmm_dem_cache_eligible(dem, ipart, jpart)`
/// checks whether the particles `ipart` and `jpart` are eligible neighbors
/// in the simulation `dem`.
/// Note that this function is neither symmetric nor reflexive
/// with respect to particle indices.
__attribute__ ((__nonnull__, __pure__))
inline bool bmm_dem_cache_eligible(struct bmm_dem const* const dem,
    size_t const ipart, size_t const jpart) {
  if (dem->cache.hurr[ipart] == dem->cache.hurr[jpart] && jpart <= ipart)
    return false;

  if (bmm_geom2d_cpdist2(dem->part.x[ipart], dem->part.x[jpart],
        dem->opts.box.x, dem->opts.box.per) >
      bmm_fp_sq(dem->opts.cache.rcutoff))
    return false;

  return true;
}

/// The call `bmm_dem_cache_x(dem)`
/// caches the position of every particle
/// in the simulation `dem`.
__attribute__ ((__nonnull__))
void bmm_dem_cache_x(struct bmm_dem*);

/// The call `bmm_dem_cache_doto(dem, ipart, mask)`
/// tries to add all the eligible particles
/// inside the `mask`-masked neighborhood
/// of the particle `ipart` to its neighbors
/// in the simulation `dem`.
/// If there is enough capacity and the operation succeeds,
/// `true` is returned.
/// Otherwise `false` is returned and
/// the previous neighbor relations are restored.
__attribute__ ((__nonnull__))
bool bmm_dem_cache_doto(struct bmm_dem*, size_t, int);

/// The call `bmm_dem_cache_dofrom(dem, ipart)`
/// tries to add the particle `ipart` to the neighbors
/// of all the eligible particles
/// inside its `mask`-masked neighborhood
/// in the simulation `dem`.
/// If there is enough capacity and the operation succeeds,
/// `true` is returned.
/// Otherwise `false` is returned and
/// the previous neighbor relations are restored.
__attribute__ ((__nonnull__))
bool bmm_dem_cache_dofrom(struct bmm_dem*, size_t, int);

/// The call `bmm_dem_ijcell(pijcell, dem, x)`
/// writes the neighbor cell of the particle at `x`
/// in the simulation `dem`
/// into the index vector `pijcell`
__attribute__ ((__nonnull__))
void bmm_dem_ijcell(size_t*, struct bmm_dem const*, double const*);

/// The call `bmm_dem_cache_cell(dem, ipart)`
/// tries to add the particle `ipart` to the appropriate neighbor cell
/// in the simulation `dem`.
/// If there is enough capacity and the operation succeeds,
/// `true` is returned.
/// Otherwise `false` is returned.
__attribute__ ((__nonnull__))
bool bmm_dem_cache_cell(struct bmm_dem*, size_t);

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

/// The call `bmm_dem_opts_def(opts)`
/// writes the default simulation options into `opts`.
__attribute__ ((__nonnull__))
void bmm_dem_opts_def(struct bmm_dem_opts*);

/// The call `bmm_dem_def(dem, opts)`
/// writes the default simulation state into `dem`
/// with the simulation options `opts`.
__attribute__ ((__nonnull__))
void bmm_dem_def(struct bmm_dem*, struct bmm_dem_opts const*);

__attribute__ ((__nonnull__))
double bmm_dem_ekinetic(struct bmm_dem const*);

__attribute__ ((__nonnull__))
double bmm_dem_pvector(struct bmm_dem const*);

__attribute__ ((__nonnull__))
double bmm_dem_pscalar(struct bmm_dem const*);

__attribute__ ((__nonnull__))
double bmm_dem_lscalar(struct bmm_dem const*);

__attribute__ ((__nonnull__))
double bmm_dem_cor(struct bmm_dem const*);

__attribute__ ((__nonnull__))
bool bmm_dem_step(struct bmm_dem*);

__attribute__ ((__nonnull__))
bool bmm_dem_comm(struct bmm_dem*);

/// The call `bmm_dem_report(dem)`
/// prints informal diagnostics for the simulation state `dem`.
__attribute__ ((__nonnull__))
bool bmm_dem_report(struct bmm_dem const*);

/// The call `bmm_dem_run(dem)`
/// processes all incoming messages and
/// handles signals with the simulation state `dem`.
__attribute__ ((__nonnull__))
bool bmm_dem_run(struct bmm_dem*);

/// The call `bmm_dem_run_with(opts)`
/// processes all incoming messages and
/// handles signals with the simulation options `opts`.
__attribute__ ((__nonnull__))
bool bmm_dem_run_with(struct bmm_dem_opts const*);

// TODO These are questionable to expose.

size_t bmm_dem_sniff_size(struct bmm_dem const* const dem,
    enum bmm_msg_num const num);

bool bmm_dem_puts_stuff(struct bmm_dem const* const dem,
    enum bmm_msg_num const num);

bool bmm_dem_puts(struct bmm_dem const* const dem,
    enum bmm_msg_num const num);

/// The call `bmm_dem_script_ongoing(dem)`
/// checks whether the simulation `dem` has not ended.
inline bool bmm_dem_script_ongoing(struct bmm_dem const* const dem) {
  return dem->script.i < dem->opts.script.n;
}

/// The call `bmm_dem_script_trans(dem)`
/// transitions the simulation `dem` to the next state and
/// checks whether the simulation did not end while doing so.
/// Make sure the simulation has not ended prior to the call
/// by calling `bmm_dem_script_ongoing`.
inline bool bmm_dem_script_trans(struct bmm_dem* const dem) {
  double const toff = dem->time.t - dem->script.tprev -
    dem->opts.script.tspan[dem->script.i];

  if (toff >= 0.0) {
    dem->script.tprev = dem->time.t;
    dem->script.ttrans[dem->script.i] = dem->time.t;
    dem->script.toff[dem->script.i] = toff;
    ++dem->script.i;

    return bmm_dem_script_ongoing(dem);
  }

  return true;
}

/// The call `bmm_dem_cache_fresh(dem)`
/// checks whether the cache for the simulation `dem` is not yet stale.
inline bool bmm_dem_cache_fresh(struct bmm_dem const* const dem) {
  for (size_t idim = 0; idim < BMM_NDIM; ++idim) {
    double const dx = dem->opts.box.x[idim] /
      (double) ((dem->opts.cache.ncell[idim] - 2) * 2);
    // TODO Use this instead of `dx - dem->part.r[ipart]`.
    // double const d = dx - dem->opts.part.rnew[1];

    for (size_t ipart = 0; ipart < dem->part.n; ++ipart)
      if (fabs(bmm_fp_swrap(dem->part.x[ipart][idim] -
              dem->cache.x[ipart][idim], dem->opts.box.x[idim])) >=
          dx - dem->part.r[ipart])
        return false;
  }

  return true;
}

#endif
