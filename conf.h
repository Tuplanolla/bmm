// Configuration constants.
#ifndef BMM_CONF_H
#define BMM_CONF_H

#include <limits.h>
#include <stddef.h>

#define BMM_MSG_MAX ((size_t) (UCHAR_MAX + 1))
#define BMM_BIN_MAX ((size_t) 1024)
#define BMM_PART_MAX ((size_t) 1024)
// #define BMM_PART_MAX ((size_t) 16384)
// Maximum number of cells per dimension.
#define BMM_CELL_MAX ((size_t) 128)
// Maximum cell or neighbor population.
#define BMM_GROUP_MAX ((size_t) 64)
#define BMM_STEP_MAX ((size_t) 16777216)

#endif
