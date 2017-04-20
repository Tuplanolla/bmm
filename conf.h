#ifndef BMM_CONF_H
/// Configuration constants.
#define BMM_CONF_H

#include <stddef.h>

// The maximum number of characters in an identifier.
#define BMM_NCHARID ((size_t) 64)

#define BMM_NMSG ((size_t) 256)
#define BMM_NBIN ((size_t) 1024)
#define BMM_NPART ((size_t) 1024)
// #define BMM_NPART ((size_t) 16384)
// Maximum number of cells per dimension.
#define BMM_NCELL ((size_t) 128)
// Maximum cell or neighbor population.
#define BMM_NGROUP ((size_t) 64)
#define BMM_NSTEP ((size_t) 16777216)

#endif
