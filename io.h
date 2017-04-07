#ifndef BMM_IO_H
/// Input and output operations.
#define BMM_IO_H

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <sys/time.h>
#include <unistd.h>

#include "ext.h"

/// This enumeration is returned by waiting operations.
enum bmm_io_wait {
  BMM_IO_WAIT_ERROR = -1,
  BMM_IO_WAIT_TIMEOUT = 0,
  BMM_IO_WAIT_READY = 1
};

/// This enumeration is returned by reading and fast-forwarding operations.
enum bmm_io_read {
  BMM_IO_READ_ERROR = -1,
  BMM_IO_READ_EOF = 0,
  BMM_IO_READ_SUCCESS = 1
};

/// The call `bmm_io_wait(fd, timeout)`
/// waits for input from the file descriptor `fd` or times out after `timeout`.
/// The remaining time is written into `timeout`.
__attribute__ ((__nonnull__))
enum bmm_io_wait bmm_io_wait(int, struct timeval*);

/// The call `bmm_io_redir(out, in, size)`
/// reads `size` bytes from `in` and writes them into `out`.
/// The return value is the number of bytes read.
__attribute__ ((__nonnull__))
size_t bmm_io_redir(FILE*, FILE*, size_t);

/// The call `bmm_io_fastfw(stream, size)`
/// reads `size` bytes from `stream` and discards them.
/// The return value is the number of bytes read.
__attribute__ ((__nonnull__))
size_t bmm_io_fastfw(FILE*, size_t);

/// The call `bmm_io_waitin(timeout)`
/// waits for input from the standard input or times out after `timeout`.
__attribute__ ((__nonnull__))
inline enum bmm_io_wait bmm_io_waitin(struct timeval* const timeout) {
  return bmm_io_wait(STDIN_FILENO, timeout);
}

/// The call `bmm_io_redirio(size)`
/// reads `size` bytes from the standard input and
/// writes them into the standard output.
inline bool bmm_io_redirio(size_t const size) {
  return bmm_io_redir(stdout, stdin, size) == size;
}

/// The call `bmm_io_fastfwin(size)`
/// reads `size` bytes from the standard input and discards them.
inline enum bmm_io_read bmm_io_fastfwin(size_t const size) {
  if (bmm_io_fastfw(stdin, size) == size)
    return BMM_IO_READ_SUCCESS;
  else if (feof(stdin) != 0)
    return BMM_IO_READ_EOF;
  else
    return BMM_IO_READ_ERROR;
}

/// The call `bmm_io_readin(ptr, size)`
/// reads `size` bytes from the standard input into `ptr`.
__attribute__ ((__nonnull__))
inline enum bmm_io_read bmm_io_readin(void* const ptr, size_t const size) {
  if (fread(ptr, size, 1, stdin) == 1)
    return BMM_IO_READ_SUCCESS;
  else if (feof(stdin) != 0)
    return BMM_IO_READ_EOF;
  else
    return BMM_IO_READ_ERROR;
}

/// The call `bmm_io_writeout(ptr, size)`
/// writes `size` bytes from `ptr` into the standard output.
__attribute__ ((__nonnull__))
inline bool bmm_io_writeout(void const* const ptr, size_t const size) {
  return fwrite(ptr, size, 1, stdout) == 1;
}

#endif
