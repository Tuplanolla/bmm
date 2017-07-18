#ifdef _GNU_SOURCE

#pragma push_macro("_GNU_SOURCE")
#undef _GNU_SOURCE

#include <string.h>

static int (*const bmm_pstrerror_r)(int, char *, size_t) = &strerror_r;

#pragma pop_macro("_GNU_SOURCE")
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#else

#include <string.h>

static int (*const bmm_pstrerror_r)(int, char *, size_t) = &strerror_r;

#endif

#include <errno.h>
#include <stdbool.h>
#include <stddef.h>

#include "hack.h"

bool bmm_hack_strerror_r(int const num, char *const buf, size_t const n) {
  int const i = (*bmm_pstrerror_r)(num, buf, n);
  if (i != 0) {
    if (i > 0)
      errno = i;

    return false;
  }

  return true;
}
