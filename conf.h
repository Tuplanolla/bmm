#ifndef BMM_CONF_H
/// Configuration constants.
#define BMM_CONF_H

#include "cpp.h"

/// Full name of the project.
#define BMM_NAME "Brittle Matter Matters"

/// Shortened name of the project.
#define BMM_SHORT "BMM"

/// Program identifier.
#define BMM_PROGID "bmm"

/// Program version numbers.
#define BMM_VERSION_MAJOR 0
#define BMM_VERSION_MINOR 0
#define BMM_VERSION_PATCH 0

/// Program version string.
#define BMM_VERSION_STRING \
  BMM_VERSION(BMM_VERSION_MAJOR, BMM_VERSION_MINOR, BMM_VERSION_PATCH)

/// Number of dimensions.
/// This is not adjustable without other changes.
#define BMM_NDIM 2

/// Maximum number of characters in an identifier.
/// This may not be adjustable without other changes.
#define BMM_NCHARID 64

/// Maximum number of particles.
#define BMM_NPART 1024

/// Maximum number of links per particle.
#define BMM_NLINK 8

/// Maximum number of neighbor cells per dimension.
#define BMM_NCELL 32

/// Maximum number of particles per neighbor cell.
#define BMM_NGROUP 64

/// Maximum number of message numbers.
/// This is not adjustable without other changes.
#define BMM_NMSG 256

/// Maximum number of simulation steps.
#define BMM_NSTEP 16777216

/// Maximum number of histogram bins.
#define BMM_NBIN 1024

#endif
