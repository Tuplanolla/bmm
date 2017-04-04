// Input and output operations.
#ifndef BMM_IO_H
#define BMM_IO_H

#include "ext.h"
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <sys/time.h>
#include <unistd.h>

enum bmm_io_wait {
  BMM_IO_WAIT_ERROR = -1,
  BMM_IO_WAIT_TIMEOUT = 0,
  BMM_IO_WAIT_READY = 1
};

__attribute__ ((__nonnull__))
enum bmm_io_wait bmm_io_wait(int, struct timeval*);

__attribute__ ((__nonnull__))
size_t bmm_io_redir(FILE*, FILE*, size_t);

__attribute__ ((__nonnull__))
size_t bmm_io_fastfw(FILE*, size_t);

__attribute__ ((__nonnull__))
inline enum bmm_io_wait bmm_io_waitin(struct timeval* const timeout) {
  return bmm_io_wait(STDIN_FILENO, timeout);
}

inline bool bmm_io_redirio(size_t const size) {
  return bmm_io_redir(stdout, stdin, size) == size;
}

enum bmm_io_read {
  BMM_IO_READ_ERROR = -1,
  BMM_IO_READ_EOF = 0,
  BMM_IO_READ_SUCCESS = 1
};

__attribute__ ((__nonnull__))
inline enum bmm_io_read bmm_io_fastfwin(size_t const size) {
  if (bmm_io_fastfw(stdin, size) == size)
    return BMM_IO_READ_SUCCESS;
  else if (feof(stdin) != 0)
    return BMM_IO_READ_EOF;
  else
    return BMM_IO_READ_ERROR;
}

__attribute__ ((__nonnull__))
inline enum bmm_io_read bmm_io_readin(void* const ptr, size_t const size) {
  if (fread(ptr, size, 1, stdin) == 1)
    return BMM_IO_READ_SUCCESS;
  else if (feof(stdin) != 0)
    return BMM_IO_READ_EOF;
  else
    return BMM_IO_READ_ERROR;
}

__attribute__ ((__nonnull__))
inline bool bmm_io_writeout(void const* const ptr, size_t const size) {
  return fwrite(ptr, size, 1, stdout) == 1;
}

#endif
