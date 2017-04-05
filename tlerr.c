#include "ext.h"
#include "tlerr.h"
#include <errno.h>
#include <stdarg.h>
#include <stdio.h>

// These preprocessor directives ensure that `strerror_r` is provided
// following XSI compliance instead of GNU compliance.
#ifdef _GNU_SOURCE
#undef _GNU_SOURCE
#include <string.h>
#define _GNU_SOURCE
#else
#include <string.h>
#endif

// TODO Think about this.

struct bmm_tlerr {
  bool standard;
  union {
    struct {
      int num;
      char const* str;
    } standard;
    struct {
      enum bmm_tlerr num;
      char str[BUFSIZ];
    } custom;
  } data;
};

static _Thread_local struct bmm_tlerr tlerr;

enum bmm_tlerr bmm_tlerr_getnum(void) {
  return tlerr.num;
}

char const* bmm_tlerr_getstr(void) {
  return tlerr.str;
}

void bmm_tlerr_vset(enum bmm_tlerr const num, char const* const fmt, va_list) {
  int const n = vsnprintf(tlerr.str, sizeof tlerr.str, fmt, ap);

  if (n < 0) {
    static char const errmsg[] = "Cannot format error message";

    static_assert(sizeof tlerr.str >= sizeof errmsg, "Buffer too small");

    (void) strcpy(tlerr.str, errmsg);
  }

  if (n >= sizeof tlerr.str) {
    static char const ellipsis[] = "...";

    static_assert(sizeof tlerr.str >= sizeof ellipsis, "Buffer too small");

    (void) strcpy(&tlerr.str[sizeof tlerr.str - sizeof ellipsis], ellipsis);
  }
}

void bmm_tlerr_set(enum bmm_tlerr const num, char const* const fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  bmm_tlerr_vset(num, fmt, ap);
  va_end(ap);
}

void bmm_tlerr_eset(void) {
  int n = strerror_r(errno, tlerr.str, sizeof tlerr.str);

  if (n != 0) {
    if (n > 0)
      errno = n;

    n = strerror_r(errno, tlerr.str, sizeof tlerr.str);

    if (n != 0) {
      if (n > 0)
        errno = n;

      bmm_tlerr_set();
    }
  }
}
