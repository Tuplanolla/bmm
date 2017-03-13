#ifndef BMM_WATCH_H
#define BMM_WATCH_H

#include "exts.h"
#include <sys/time.h>

enum watch {
  WATCH_ERROR,
  WATCH_TIMEOUT,
  WATCH_READY
};

__attribute__ ((__nonnull__))
enum watch watchin(struct timeval*);

#endif
