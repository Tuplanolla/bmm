/// Sophisticated error reporting.
#ifndef BMM_ERR_H
#define BMM_ERR_H

#include <stdarg.h>
#include <stddef.h>
#include <stdlib.h>

#include "ext.h"

/// The call `bmm_err_reset(str)` initializes the timing mechanism
/// `BMM_ERR_WARN` and `BMM_ERR_ABORT` are based on for a program named `str`.
/// Calling `bmm_err_reset` before using any other error handling procedures
/// from this compilation unit is not necessary.
/// If `bmm_err_reset` is never called,
/// the reference point for timestamps in error messages is chosen arbitrarily.
void bmm_err_reset(char const*);

/// The call `bmm_err_warn(func, file, line, str)`
/// prints a detailed error message
/// from the procedure `func` on line `line` of file `file`
/// to the standard error output stream,
/// where the pointer `str` may be `NULL` or point to a string
/// that contains the name of the procedure or variable that caused the error.
/// The error message is taken from the external variable `errno`,
/// which may be reset in the process.
///
/// The format of the error message consists of
///
/// * `"[time]"`, mimicking `dmesg`,
/// * `"file:line"`, mimicking `cc` and
/// * `"procedure: message"`, mimicking `perror`.
///
/// An example message follows (sans wrapping).
///
///     [1.048576]
///     bmm-dem (16384):
///     bmm_dem_run (dem.c:420):
///     malloc: Cannot allocate memory
///
/// Note that it may be more convenient to use `BMM_ERR_WARN`
/// instead of calling `bmm_err_warn` directly.
__attribute__ ((__nonnull__ (1, 2)))
void bmm_err_warn(char const*, char const*, size_t, char const*);

/// The call `bmm_err_vfwarn(func, file, line, str, fmt, ap)`
/// prints a detailed error message
/// from the procedure `func` on line `line` of file `file`
/// to the standard error output stream,
/// where the pointer `str` may be `NULL` or point to a string
/// that contains the name of the procedure or variable that caused the error.
/// The error message is built with the format string `fmt` and arguments `ap`.
///
/// See `bmm_err_warn`.
__attribute__ ((__format__ (__printf__, 5, 0), __nonnull__ (1, 2, 5)))
void bmm_err_vfwarn(char const*, char const*,
    size_t, char const*, char const*, va_list);

/// See `bmm_err_vfwarn`.
__attribute__ ((__format__ (__printf__, 5, 6), __nonnull__ (1, 2, 5)))
void bmm_err_fwarn(char const*, char const*,
    size_t, char const*, char const*, ...);

/// The preprocessor directive `BMM_ERR_WARN(ptr)`
/// prints a detailed error message
/// to the standard error output stream,
/// where the pointer `ptr` may be `NULL` or
/// point to the procedure or variable that caused the error.
///
/// See `bmm_err_warn(ptr)`.
#define BMM_ERR_WARN(ptr) (bmm_err_warn(__func__, __FILE__, \
      (size_t) __LINE__, (ptr) == NULL ? NULL : #ptr))

/// See `BMM_ERR_WARN(ptr)`.
#define BMM_ERR_VFWARN(ptr, fmt, ap) (bmm_err_vfwarn(__func__, __FILE__, \
      (size_t) __LINE__, (ptr) == NULL ? NULL : #ptr, fmt, ap))

/// See `BMM_ERR_VFWARN(ptr)`.
#define BMM_ERR_FWARN(ptr, ...) (bmm_err_fwarn(__func__, __FILE__, \
      (size_t) __LINE__, (ptr) == NULL ? NULL : #ptr, __VA_ARGS__))

/// The call `BMM_ERR_ABORT(ptr)` is equivalent
/// to `BMM_ERR_WARN(ptr)` followed by `abort()`.
#define BMM_ERR_ABORT(ptr) (BMM_ERR_WARN(ptr), abort())

/// See `BMM_ERR_ABORT`.
#define BMM_ERR_VFABORT(ptr, fmt, ap) (BMM_ERR_VFWARN(ptr, fmt, ap), abort())

/// See `BMM_ERR_VFABORT`.
#define BMM_ERR_FABORT(ptr, ...) (BMM_ERR_FWARN(ptr, __VA_ARGS__), abort())

#endif
