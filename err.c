#include "err.h"
#include "sec.h"
#include <errno.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <unistd.h>

static double t0 = 0.0;

void bmm_err_reset(void) {
  t0 = bmm_sec_now();
}

void bmm_err_warn(char const* const func, char const* const file,
    size_t const line, char const* const str) {
  int const tmp = errno;

  double const t1 = bmm_sec_now();

  if (fprintf(stderr, "[%f] (%zu) <%s> %s:%zu: ",
      t1 - t0, (size_t) getpid(), func, file, line) < 0)
    return;

  errno = tmp;

  perror(str);
}

void bmm_err_vfwarn(char const* const func, char const* const file,
    size_t const line, char const* const str,
    char const* const fmt, va_list ap) {
  double const t1 = bmm_sec_now();

  if (fprintf(stderr, "[%f] (%zu) <%s> %s:%zu: ",
      t1 - t0, (size_t) getpid(), func, file, line) < 0)
    return;

  if (str != NULL)
    if (fprintf(stderr, "%s: ", str) < 0)
      return;

  if (vfprintf(stderr, fmt, ap) < 0)
    return;

  (void) fputc('\n', stderr);
}

void bmm_err_fwarn(char const* const func, char const* const file,
    size_t const line, char const* const str,
    char const* const fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  bmm_err_vfwarn(func, file, line, str, fmt, ap);
  va_end(ap);
}
