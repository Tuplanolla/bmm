#include <errno.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "err.h"
#include "sec.h"

static char const* prog = NULL;

static double t0 = 0.0;

void bmm_err_reset(char const* const str) {
  char const* const chr = strrchr(str, '/');
  prog = chr == NULL ? prog : &chr[1];

  t0 = bmm_sec_now();
}

void bmm_err_warn(char const* const func, char const* const file,
    size_t const line, char const* const str) {
  int const tmp = errno;

  if (prog == NULL)
#ifdef _GNU_SOURCE
    prog = program_invocation_name;
#else
    prog = "a.out";
#endif

  double const t1 = bmm_sec_now();

  if (fprintf(stderr, "[%f] %s (%zu): %s (%s:%zu): ",
      t1 - t0, prog, (size_t) getpid(), func, file, line) < 0)
    return;

  errno = tmp;

  perror(str);
}

void bmm_err_vfwarn(char const* const func, char const* const file,
    size_t const line, char const* const str,
    char const* const fmt, va_list ap) {
  double const t1 = bmm_sec_now();

  if (prog == NULL)
#ifdef _GNU_SOURCE
    prog = program_invocation_name;
#else
    prog = "a.out";
#endif

  if (fprintf(stderr, "[%f] %s (%zu): %s (%s:%zu): ",
      t1 - t0, prog, (size_t) getpid(), func, file, line) < 0)
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
