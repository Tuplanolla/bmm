// Support for thread-local error handling.
#ifndef BMM_TLERR_H
#define BMM_TLERR_H

// TODO Think about this.
// TODO Get rid of conflicting `err.h`.

#include "ext.h"
#include <stdarg.h>

enum bmm_tlerr {
  BMM_TLERR_SUCCESS = 0,
  BMM_TLERR_ERRMSG = 1,
  BMM_TLERR_BAD,
  BMM_TLERR_VERYBAD
};

enum bmm_tlerr bmm_tlerr_getnum(void);

char const* bmm_tlerr_getstr(void);

__attribute__ ((__format__ (__printf__, 2, 0), __nonnull__))
void bmm_tlerr_vset(enum bmm_tlerr, char const*, va_list);

__attribute__ ((__format__ (__printf__, 2, 3), __nonnull__))
void bmm_tlerr_set(enum bmm_tlerr, char const*, ...);

void bmm_tlerr_eset(void);
