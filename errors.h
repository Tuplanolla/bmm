#ifndef BMM_ERRORS_H
#define BMM_ERRORS_H

#include "exts.h"
#include <stdarg.h>

// The call `bmm_verror(fmt, ap)` prints an error message
// with the format `fmt` and arguments `ap`.
__attribute__ ((__format__ (__printf__, 1, 0), __nonnull__))
void bmm_verror(char const*, va_list);

// The call `bmm_error(fmt, x, y, z)` prints an error message
// with the format `fmt` and arguments `x`, `y` and `z`.
// Any other number of arguments works too.
__attribute__ ((__format__ (__printf__, 1, 2), __nonnull__))
void bmm_error(char const*, ...);

#endif
