// String operations.
#ifndef BMM_STR_H
#define BMM_STR_H

#include <stdbool.h>
#include <stddef.h>

#include "ext.h"

// The call `bmm_str_strtou(ptr, str)`
// parses an unsigned integer from the string `str` and saves it into `ptr`.
__attribute__ ((__nonnull__ (2)))
bool bmm_str_strtou(unsigned int*, char const*);

// The call `bmm_str_strtoz(ptr, str)`
// parses a size from the string `str` and saves it into `ptr`.
__attribute__ ((__nonnull__ (2)))
bool bmm_str_strtoz(size_t*, char const*);

// The call `bmm_str_strtod(ptr, str)`
// parses a floating-point number from the string `str` and
// saves it into `ptr`.
__attribute__ ((__nonnull__ (2)))
bool bmm_str_strtod(double*, char const*);

// The call `bmm_str_strtob(ptr, str)`
// parses a truth value from the string `str` and
// saves it into `ptr`.
__attribute__ ((__nonnull__ (2)))
bool bmm_str_strtob(bool*, char const*);

#endif
