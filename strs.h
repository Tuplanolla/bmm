#ifndef BMM_STRS_H
#define BMM_STRS_H

#include "exts.h"
#include <stdbool.h>
#include <stddef.h>

// The call `bmm_strtoz(str, p)`
// parses a size from the string `str` and saves it into `p`.
__attribute__ ((__nonnull__ (1)))
bool bmm_strtoz(char const*, size_t*);

#endif
