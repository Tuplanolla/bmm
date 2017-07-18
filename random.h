#ifndef BMM_RANDOM_H
/// Random number generation.
#define BMM_RANDOM_H

#include <gsl/gsl_rng.h>

#include "ext.h"

__attribute__ ((__nonnull__))
inline double bmm_random_get(gsl_rng *const rng, double const *const x) {
  return gsl_rng_uniform(rng) * (x[1] - x[0]) + x[0];
}

#endif
