// String operations.
#ifndef BMM_STR_H
#define BMM_STR_H

#include "ext.h"
#include <stdbool.h>
#include <stddef.h>

// The call `bmm_str_strtoz(str, ptr)`
// parses a size from the string `str` and saves it into `ptr`.
__attribute__ ((__nonnull__ (1)))
bool bmm_str_strtoz(char const*, size_t*);

#endif
