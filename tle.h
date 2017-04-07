// Support for thread-local error handling.
#ifndef BMM_TLE_H
#define BMM_TLE_H

#include <stdarg.h>
#include <stddef.h>

#include "ext.h"

// This enumeration complements std thread-local error numbers.
enum bmm_tle {
  BMM_TLE_SUCCESS = 0,
  BMM_TLE_BAD = 1,
  BMM_TLE_VERYBAD
};

// The call `bmm_tle_reset(str)`
// initializes the error handling mechanism for a program named `str`.
// Calling `bmm_tle_reset` before using any of the error handling procedures
// from this compilation unit is not necessary,
// but neglecting to do so will result in a default program name and
// an arbitrary timestamp offset.
//
// Note that as of now this is only compatible with ASCII strings.
void bmm_tle_reset(char const*);

// The call `bmm_tle_num_std()`
// returns the standard thread-local error number
// if it was set with `bmm_tle_std`.
// Otherwise success is returned.
__attribute__ ((__pure__))
int bmm_tle_num_std(void);

// The call `bmm_tle_num_ext()`
// returns the extended thread-local error number
// if it was set with `bmm_tle_ext` or `bmm_tle_vext`.
// Otherwise success is returned.
__attribute__ ((__pure__))
enum bmm_tle bmm_tle_num_ext(void);

// The call `bmm_tle_msg()`
// returns the thread-local error message.
__attribute__ ((__pure__))
char const* bmm_tle_msg(void);

// The call `bmm_tle_std()`
// sets the thread-local error number to the standard error number `errno` and
// copies the corresponding error message
// with the assistance of `strerror_r`.
void bmm_tle_std(void);

// The call `bmm_tle_stds(file, line, proc)`
// sets the thread-local error number to the standard error number `errno` and
// copies the corresponding error message
// for the procedure `proc` on line `line` of file `file`
// with the assistance of `strerror_r`.
//
// The error message is formatted with
//
// * `"[time]"`, mimicking `dmesg`,
// * `"file:line"`, mimicking `cc` and
// * `"procedure: message"`, mimicking `perror`
//
// and looks as follows (except for wrapping).
//
//     [1.048576]
//     bmm-dem (16384):
//     bmm_dem_run (dem.c:420):
//     malloc: Cannot allocate memory
//
// Note that it may be more convenient to use `BMM_TLE_STDS`
// instead of calling `bmm_tle_stds` directly.
__attribute__ ((__nonnull__))
void bmm_tle_stds(char const*, size_t, char const*);

// See `bmm_tle_stds`.
#define BMM_TLE_STDS() \
  (bmm_tle_stds(__FILE__, (size_t) __LINE__, __func__))

// The call `bmm_tle_vext(num, fmt, ap)`
// sets the thread-local error number to the extended error number `num` and
// builds the corresponding error message
// with the format string `fmt` and arguments `ap`.
__attribute__ ((__format__ (__printf__, 2, 0), __nonnull__))
void bmm_tle_vext(enum bmm_tle, char const*, va_list);

// See `bmm_tle_vext`.
__attribute__ ((__format__ (__printf__, 2, 3), __nonnull__))
void bmm_tle_ext(enum bmm_tle, char const*, ...);

// The call `bmm_tle_vext(num, file, line, proc, fmt, ap)`
// sets the thread-local error number to the extended error number `num` and
// builds the corresponding error message
// for the procedure `proc` on line `line` of file `file`
// with the format string `fmt` and arguments `ap`.
//
// See `bmm_tle_stds`.
__attribute__ ((__format__ (__printf__, 5, 0), __nonnull__))
void bmm_tle_vexts(char const*, size_t, char const*,
    enum bmm_tle, char const*, va_list);

// See `bmm_tle_vexts`.
#define BMM_TLE_VEXTS(num, fmt, ap) \
  (bmm_tle_vexts(__FILE__, (size_t) __LINE__, __func__, num, fmt, ap))

// See `bmm_tle_vexts`.
__attribute__ ((__format__ (__printf__, 5, 6), __nonnull__))
void bmm_tle_exts(char const*, size_t, char const*,
    enum bmm_tle, char const*, ...);

// See `bmm_tle_exts`.
#define BMM_TLE_EXTS(num, fmt, ...) \
  (bmm_tle_exts(__FILE__, (size_t) __LINE__, __func__, num, fmt, __VA_ARGS__))

#endif
