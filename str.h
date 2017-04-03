// String operations.
#ifndef BMM_STR_H
#define BMM_STR_H

#include "ext.h"
#include <stdbool.h>
#include <stddef.h>

// The call `bmm_str_strtou(ptr, str)`
// parses an unsigned integer from the string `str` and saves it into `ptr`.
__attribute__ ((__nonnull__ (1)))
bool bmm_str_strtou(unsigned int*, char const*);

// The call `bmm_str_strtoz(ptr, str)`
// parses a size from the string `str` and saves it into `ptr`.
__attribute__ ((__nonnull__ (1)))
bool bmm_str_strtoz(size_t*, char const*);

// The call `bmm_str_strtod(ptr, str)`
// parses a floating-point number from the string `str` and
// saves it into `ptr`.
__attribute__ ((__nonnull__ (1)))
bool bmm_str_strtod(double*, char const*);

#endif
