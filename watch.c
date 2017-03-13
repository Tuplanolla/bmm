#include "watch.h"
#include <sys/select.h>
#include <unistd.h>

enum watch watchin(struct timeval* const timeout) {
  fd_set fds;
  FD_ZERO(&fds);
  FD_SET(STDIN_FILENO, &fds);

  int const result = select(1, &fds, NULL, NULL, timeout);
  if (result == -1)
    return WATCH_ERROR;
  else if (result == 0)
    return WATCH_TIMEOUT;
  else
    return WATCH_READY;
}

