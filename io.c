#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/select.h>

#include "io.h"
#include "tle.h"

extern inline bool bmm_io_read_to_bool(enum bmm_io_read);

enum bmm_io_wait bmm_io_wait(int const fd, struct timeval* const timeout) {
  fd_set fds;
  FD_ZERO(&fds);
  FD_SET(fd, &fds);

  int const result = select(1, &fds, NULL, NULL, timeout);
  switch (result) {
    case -1:
      return BMM_IO_WAIT_ERROR;
    case 0:
      return BMM_IO_WAIT_TIMEOUT;
    default:
      return BMM_IO_WAIT_READY;
  }
}

size_t bmm_io_redir(FILE* const out, FILE* const in, size_t const size) {
  size_t progress = 0;

  unsigned char buf[BUFSIZ];

  while (progress < size) {
    size_t const ndiff = size - progress;
    size_t const nmemb = ndiff < sizeof buf ? ndiff : sizeof buf;

    size_t const nread = fread(buf, 1, nmemb, in);
    if (nread == 0)
      break;

    size_t const nwritten = fwrite(buf, 1, nmemb, out);
    if (nwritten < nread)
      break;

    progress += nread;
  }

  return progress;
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

enum bmm_io_wait bmm_io_waitin(struct timeval*);

extern inline bool bmm_io_redirio(size_t);

extern inline enum bmm_io_read bmm_io_fastfwin(size_t);

extern inline enum bmm_io_read bmm_io_readin(void*, size_t);

extern inline bool bmm_io_writeout(void const*, size_t);

void* bmm_io_fconts(size_t* const ptr, FILE* const stream) {
  if (fseek(stream, 0, SEEK_END) != 0) {
    BMM_TLE_STDS();

    return NULL;
  }

  long int const pos = ftell(stream);

  if (fseek(stream, 0, SEEK_SET) != 0) {
    BMM_TLE_STDS();

    return NULL;
  }

  size_t const size = (size_t) pos;

  void* const buf = malloc(size);
  if (buf == NULL) {
    BMM_TLE_STDS();

    return NULL;
  }

  if (fread(buf, size, 1, stream) != 1) {
    BMM_TLE_STDS();

    free(buf);

    return NULL;
  }

  *ptr = size;

  return buf;
}

void* bmm_io_conts(size_t* const ptr, char const* const path) {
  FILE* const stream = fopen(path, "r");
  if (stream == NULL) {
    BMM_TLE_STDS();

    return NULL;
  }

  void* const buf = bmm_io_fconts(ptr, stream);
  if (buf == NULL) {
    BMM_TLE_STDS();

    if (fclose(stream) != 0)
      BMM_TLE_STDS();

    return NULL;
  }

  if (fclose(stream) != 0) {
    BMM_TLE_STDS();

    free(buf);

    return NULL;
  }

  return buf;
}
