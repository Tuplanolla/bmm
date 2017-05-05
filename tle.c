#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "ext.h"
#include "hack.h"
#include "sec.h"
#include "tle.h"

struct bmm_tle {
  char const* prog;
  double sec;
  enum bmm_tle_tag tag;
  union {
    int std;
    enum bmm_tle_num ext;
  } num;
  char buf[BUFSIZ];
};

static thread_local struct bmm_tle tle = {
  .prog = NULL,
  .sec = (double) NAN,
  .tag = BMM_TLE_TAG_EXT,
  .num.std = BMM_TLE_NUM_SUCCESS,
  .buf = "Success"
};

void bmm_tle_reset(char const* const prog) {
  char const* const str = strrchr(prog, '/');

  tle.prog = str == NULL ? prog : &str[1];

  tle.sec = bmm_sec_now();
}

char const* bmm_tle_prog(void) {
  return tle.prog;
}

enum bmm_tle_tag bmm_tle_tag(void) {
  return tle.tag;
}

int bmm_tle_num_std(void) {
  switch (tle.tag) {
    case BMM_TLE_TAG_STD:
      return tle.num.std;
    case BMM_TLE_TAG_EXT:
      return 0;
  }

  dynamic_assert(false, "Nonexhaustive cases");
}

enum bmm_tle_num bmm_tle_num_ext(void) {
  switch (tle.tag) {
    case BMM_TLE_TAG_STD:
      return BMM_TLE_NUM_SUCCESS;
    case BMM_TLE_TAG_EXT:
      return tle.num.ext;
  }

  dynamic_assert(false, "Nonexhaustive cases");
}

char const* bmm_tle_msg(void) {
  return tle.buf;
}

void bmm_tle_fput(FILE* const stream) {
  (void) fprintf(stream, "%s\n", tle.buf);
}

void bmm_tle_put(void) {
  (void) fprintf(stderr, "%s\n", tle.buf);
}

static void init(void) {
  if (tle.prog == NULL) {
    tle.prog = "a.out";

    tle.sec = bmm_sec_now();
  }
}

__attribute__ ((__nonnull__ (2, 4)))
static bool prefix(size_t* const ptr,
    char const* const file, size_t const line, char const* const proc) {
  init();

  int const i = snprintf(tle.buf, sizeof tle.buf,
      "[%f] %s (%zu): %s (%s:%zu): ",
      bmm_sec_now() - tle.sec, tle.prog, (size_t) getpid(), proc, file, line);

  if (i < 0)
    return false;

  size_t const n = (size_t) i;

  if (n >= sizeof tle.buf)
    return false;

  if (ptr != NULL)
    *ptr = n;

  return true;
}

static void suffix(void) {
  static char const buf[] = "Cannot report error";
  static_assert(sizeof tle.buf >= sizeof buf, "Buffer too short");
  (void) strcpy(tle.buf, buf);
}

static bool suffix_std(size_t const n) {
  tle.tag = BMM_TLE_TAG_STD;

  tle.num.std = errno;

  if (!bmm_hack_strerror_r(errno, &tle.buf[n], sizeof tle.buf - n)) {
    tle.num.std = errno;

    if (!bmm_hack_strerror_r(errno, &tle.buf[n], sizeof tle.buf - n)) {
      tle.num.std = errno;

      return false;
    }
  }

  return true;
}

__attribute__ ((__format__ (__printf__, 3, 0), __nonnull__))
static bool suffix_ext(size_t const n, enum bmm_tle_num const num,
    char const* const fmt, va_list ap) {
  tle.tag = BMM_TLE_TAG_EXT;

  tle.num.ext = num;

  size_t const m = sizeof tle.buf - n;

  int const i = vsnprintf(&tle.buf[n], m, fmt, ap);

  if (i < 0)
    return false;

  size_t const k = (size_t) i;

  if (k >= m)
    return false;

  return true;
}

void bmm_tle_std(void) {
  if (!suffix_std(0))
    suffix();
}

void bmm_tle_stds(char const* const file, size_t const line,
    char const* const proc) {
  size_t n;

  if (!prefix(&n, file, line, proc)) {
    bmm_tle_std();

    return;
  }

  if (!suffix_std(n))
    suffix();
}

void bmm_tle_vext(enum bmm_tle_num const num,
    char const* const fmt, va_list ap) {
  tle.tag = BMM_TLE_TAG_EXT;

  tle.num.ext = num;

  int const i = vsnprintf(tle.buf, sizeof tle.buf, fmt, ap);
  if (i < 0 || (size_t) i >= sizeof tle.buf)
    bmm_tle_std();
}

void bmm_tle_ext(enum bmm_tle_num const num, char const* const fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  bmm_tle_vext(num, fmt, ap);
  va_end(ap);
}

void bmm_tle_vexts(char const* const file, size_t const line,
    char const* const proc, enum bmm_tle_num const num,
    char const* const fmt, va_list ap) {
  size_t n;

  if (!prefix(&n, file, line, proc)) {
    if (!suffix_std(0))
      suffix();

    return;
  }

  if (!suffix_ext(n, num, fmt, ap))
    if (!suffix_std(n))
      suffix();
}

void bmm_tle_exts(char const* const file, size_t const line,
    char const* const proc, enum bmm_tle_num const num,
    char const* const fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  bmm_tle_vexts(file, line, proc, num, fmt, ap);
  va_end(ap);
}
