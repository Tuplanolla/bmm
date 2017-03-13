// Input and output operations.
#ifndef BMM_IO_H
#define BMM_IO_H

#include "ext.h"
#include <sys/time.h>

enum bmm_io {
  BMM_IO_ERROR = -1,
  BMM_IO_TIMEOUT = 0,
  BMM_IO_READY = 1
};

__attribute__ ((__nonnull__))
enum bmm_io bmm_io_waitin(struct timeval*);

#endif
