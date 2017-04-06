#include <errno.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "ext.h"
#include "hack.h"
#include "sec.h"
#include "tle.h"

struct tle {
  char const* prog;
  double sec;
  enum {STD, EXT} tag;
  union {
    int std;
    enum bmm_tle ext;
  } num;
  char buf[70];
  // char buf[BUFSIZ];
};

// TODO Reduce repetition and horrible control flow.
// TODO Rethink order afterwards.

static _Thread_local struct tle tle = {
  .prog = NULL,
  .sec = 0.0,
  .tag = STD,
  .num.std = 0
};

static void panic(void) {
  static char const str[] = "Cannot report error";

  static_assert(sizeof tle.buf >= sizeof str, "Buffer too short");

  (void) strcpy(tle.buf, str);
}

static void prepare(void) {
  if (tle.prog == NULL) {
    tle.prog = "a.out";

    tle.sec = bmm_sec_now();
  }
}

void bmm_tle_reset(char const* const prog) {
  char const* const str = strrchr(prog, '/');

  tle.prog = str == NULL ? prog : &str[1];

  tle.sec = bmm_sec_now();
}

int bmm_tle_num_std(void) {
  switch (tle.tag) {
    case STD:
      return tle.num.std;
    case EXT:
      return 0;
  }
}

enum bmm_tle bmm_tle_num_ext(void) {
  switch (tle.tag) {
    case STD:
      return BMM_TLE_SUCCESS;
    case EXT:
      return tle.num.ext;
  }
}

char const* bmm_tle_msg(void) {
  return tle.buf;
}

static int common_s(char const* const file, size_t const line,
    char const* const proc) {
  prepare();

  double const sec = bmm_sec_now();

  return snprintf(tle.buf, sizeof tle.buf,
      "[%f] %s (%zu): %s (%s:%zu): ",
      sec - tle.sec, tle.prog, (size_t) getpid(), proc, file, line);
}

// The basic idea for long-formatting:
// try to long-format prefix
// on error try to short-format sprintf error
//          on error try to short-format strerror error
//                   on error short-format strerror error
// try to short-format suffix
// on error try to short-format suffix sprintf error
//          on error try to short-format suffix strerror error
//                   on error short-format strerror error.

void bmm_tle_stds(char const* const file, size_t const line,
    char const* const proc) {
  tle.tag = STD;

  int const j = common_s(file, line, proc);
  if (j < 0) {
    bmm_tle_std();

    return;
  }

  size_t const n = (size_t) j;

  if (n >= sizeof tle.buf) {
    bmm_tle_std();

    return;
  }

  size_t const m = sizeof tle.buf - n;

  // Diff point.

  size_t tries = 2;

try:
  tle.num.std = errno;

  if (tries == 0)
    panic();
  else {
    if (!bmm_hack_strerror_r(errno, &tle.buf[n], m)) {
      --tries;

      goto try;
    }
  }
}

void bmm_tle_std(void) {
  tle.tag = STD;

  size_t tries = 2;

try:
  tle.num.std = errno;

  if (tries == 0)
    panic();
  else {
    if (!bmm_hack_strerror_r(errno, tle.buf, sizeof tle.buf)) {
      --tries;

      goto try;
    }
  }
}

void bmm_tle_vexts(char const* const file, size_t const line,
    char const* const proc,
    enum bmm_tle const num, char const* const fmt, va_list ap) {
  tle.tag = EXT;

  int const i = common_s(file, line, proc);
  if (i < 0) {
    bmm_tle_std();

    return;
  }

  size_t const n = (size_t) i;

  if (n >= sizeof tle.buf) {
    bmm_tle_std();

    return;
  }

  size_t const m = sizeof tle.buf - n;

  int const j = vsnprintf(&tle.buf[n], m, fmt, ap);
  if (j < 0) {
    BMM_TLE_STDS();

    return;
  }

  size_t const k = (size_t) j;

  if (k >= m) {
    BMM_TLE_STDS();

    return;
  }

  // Diff point.

  tle.num.ext = num;
}

void bmm_tle_exts(char const* const file, size_t const line,
    char const* const proc,
    enum bmm_tle const num, char const* const fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  bmm_tle_vexts(file, line, proc, num, fmt, ap);
  va_end(ap);
}

void bmm_tle_vext(enum bmm_tle const num, char const* const fmt, va_list ap) {
  tle.tag = EXT;

  tle.num.ext = num;

  int const i = vsnprintf(tle.buf, sizeof tle.buf, fmt, ap);
  if (i < 0 || (size_t) i >= sizeof tle.buf)
    bmm_tle_std();
}

void bmm_tle_ext(enum bmm_tle const num, char const* const fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  bmm_tle_vext(num, fmt, ap);
  va_end(ap);
}
