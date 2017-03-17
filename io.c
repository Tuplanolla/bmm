#include "io.h"
#include <stddef.h>
#include <stdio.h>
#include <sys/select.h>
#include <unistd.h>

enum bmm_io bmm_io_waitin(struct timeval* const timeout) {
  fd_set fds;
  FD_ZERO(&fds);
  FD_SET(STDIN_FILENO, &fds);

  int const result = select(1, &fds, NULL, NULL, timeout);
  switch (result) {
    case -1:
      return BMM_IO_ERROR;
    case 0:
      return BMM_IO_TIMEOUT;
    default:
      return BMM_IO_READY;
  }
}

size_t bmm_io_fastfw(FILE* const stream, size_t const size) {
  size_t progress = 0;

  unsigned char buf[BUFSIZ];

  while (progress < size) {
    size_t const ndiff = size - progress;
    size_t const nmemb = ndiff < sizeof buf ? ndiff : sizeof buf;

    size_t const nread = fread(buf, 1, nmemb, stream);
    if (nread == 0)
      break;

    progress += nread;
  }

  return progress;
}
