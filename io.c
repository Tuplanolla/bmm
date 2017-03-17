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

size_t bmm_io_fastfw(size_t const size, size_t const nmemb,
    FILE* const stream) {
  size_t nsum = 0;

  unsigned char buf[BUFSIZ];

  while (nsum < nmemb) {
    size_t const ndiff = nmemb - nsum;
    size_t const nbuf = ndiff < sizeof buf ? ndiff : sizeof buf;

    size_t const nread = fread(buf, size, nbuf, stream);
    if (nread == 0)
      break;

    nsum += nread;
  }

  return nsum;
}
