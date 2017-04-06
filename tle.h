// Support for thread-local error handling.
#ifndef BMM_TLE_H
#define BMM_TLE_H

#include <stdarg.h>

#include "ext.h"

// This enumeration complements std thread-local error numbers.
enum bmm_tle {
  BMM_TLE_SUCCESS = 0,
  BMM_TLE_BAD = 1,
  BMM_TLE_VERYBAD
};

// The call `bmm_tle_num_std()`
// returns the thread-local error number
// if it was set with `bmm_tle_std`.
// Otherwise success is returned.
__attribute__ ((__pure__))
int bmm_tle_num_std(void);

// The call `bmm_tle_num_ext()`
// returns the thread-local error number
// if it was set with `bmm_tle_ext` or `bmm_tle_vext`.
// Otherwise success is returned.
__attribute__ ((__pure__))
enum bmm_tle bmm_tle_num_ext(void);

// The call `bmm_tle_msg()`
// returns the thread-local error message.
__attribute__ ((__pure__))
char const* bmm_tle_msg(void);

// The call `bmm_tle_vext(num, fmt, ap)`
// sets the thread-local error number to the extended error number `num` and
// builds the corresponding error message
// with the format string `fmt` and arguments `ap`.
__attribute__ ((__format__ (__printf__, 2, 0), __nonnull__))
void bmm_tle_vext(enum bmm_tle, char const*, va_list);

// See `bmm_tle_vext`.
__attribute__ ((__format__ (__printf__, 2, 3), __nonnull__))
void bmm_tle_ext(enum bmm_tle, char const*, ...);

// The call `bmm_tle_std()`
// sets the thread-local error number to the standard error number `errno` and
// copies the corresponding error message with `strerror_r`.
void bmm_tle_std(void);

#endif
