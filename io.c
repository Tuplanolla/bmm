#include "io.h"
#include <sys/select.h>
#include <unistd.h>

enum bmm_io bmm_io_waitin(struct timeval* const timeout) {
  fd_set fds;
  FD_ZERO(&fds);
  FD_SET(STDIN_FILENO, &fds);

  int const result = select(1, &fds, NULL, NULL, timeout);
  if (result == -1)
    return BMM_IO_ERROR;
  else if (result == 0)
    return BMM_IO_TIMEOUT;
  else
    return BMM_IO_READY;
}

