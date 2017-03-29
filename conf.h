// Configuration constants.
#ifndef BMM_CONF_H
#define BMM_CONF_H

#include <stddef.h>

#define BMM_BIN_MAX ((size_t) 1024)
#define BMM_PART_MAX ((size_t) 256)
// #define BMM_PART_MAX ((size_t) 16384)
// Maximum number of cells per dimension.
#define BMM_CELL_MAX ((size_t) 64)
// Maximum cell or neighbor population.
#define BMM_GROUP_MAX ((size_t) 32)
#define BMM_STEP_MAX ((size_t) 16777216)

#endif
