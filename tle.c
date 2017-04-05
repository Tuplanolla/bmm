// These preprocessor directives try to get `strerror_r`
// with XSI compliance instead of GNU compliance.
// They have to be the very first thing in this file.
#ifdef _GNU_SOURCE
#define BMM_TLE_GNU_SOURCE _GNU_SOURCE
#undef _GNU_SOURCE
#include <string.h>
#define _GNU_SOURCE BMM_TLE_GNU_SOURCE
#undef BMM_TLE_GNU_SOURCE
#else
#include <string.h>
#endif

// This redundant prototype ensures that the attempt was successful.
int strerror_r(int, char*, size_t);

#include "ext.h"
#include "tle.h"
#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

enum tag {
  STANDARD,
  EXTENDED
};

struct tle {
  enum tag tag;
  union {
    int std;
    enum bmm_tle ext;
  } num;
  char buf[BUFSIZ];
};

static _Thread_local struct tle tle;

static void panic(void) {
  static char const str[] = "Cannot report error";

  static_assert(sizeof tle.buf >= sizeof str, "Buffer too short");

  (void) strcpy(tle.buf, str);
}

int bmm_tle_num_std(void) {
  switch (tle.tag) {
    case STANDARD:
      return tle.num.std;
    case EXTENDED:
      return 0;
  }
}

enum bmm_tle bmm_tle_num_ext(void) {
  switch (tle.tag) {
    case STANDARD:
      return BMM_TLE_SUCCESS;
    case EXTENDED:
      return tle.num.ext;
  }
}

char const* bmm_tle_msg(void) {
  return tle.buf;
}

void bmm_tle_std(void) {
  tle.tag = STANDARD;

  size_t tries = 2;

try:
  tle.num.std = errno;

  if (tries == 0)
    panic();
  else {
    int const n = strerror_r(errno, tle.buf, sizeof tle.buf);
    if (n != 0) {
      if (n > 0)
        errno = n;

      --tries;

      goto try;
    }
  }
}

void bmm_tle_vext(enum bmm_tle const num, char const* const fmt, va_list ap) {
  tle.tag = EXTENDED;

  tle.num.ext = num;

  int const n = vsnprintf(tle.buf, sizeof tle.buf, fmt, ap);
  if (n < 0 || n >= sizeof tle.buf)
    bmm_tle_std();
}

void bmm_tle_ext(enum bmm_tle const num, char const* const fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  bmm_tle_vext(num, fmt, ap);
  va_end(ap);
}
