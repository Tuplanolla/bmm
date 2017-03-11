#include "errors.h"
#include <stdarg.h>
#include <stdio.h>

void bmm_verror(char const* const fmt, va_list ap) {
  (void) vfprintf(stderr, fmt, ap);
  (void) fputc('\n', stderr);
}

void bmm_error(char const* const fmt, ...) {
  va_list ap;
  va_start(ap, fmt);

  bmm_verror(fmt, ap);

  va_end(ap);
}
