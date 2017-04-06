#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#include "bit.h"
#include "cpp.h"
#include "dem.h"
#include "err.h"
#include "io.h"
#include "msg.h"
#include "size.h"

void bmm_head_def(struct bmm_msg_head* const head) {
  head->flags = 0;
  head->type = 0;
}

// TODO Put this elsewhere.
bool bmm_msg_precheck(size_t const size, size_t const bodysize) {
  if (size < bodysize) {
    BMM_ERR_FWARN(bmm_size_from_buffer, "Message body would overflow");

    return false;
  }

  if (size != bodysize)
    BMM_ERR_FWARN(bmm_size_from_buffer, "Message body size mismatch");

  return true;
}

// Assuming `head` has already been read...
// TODO Passing `size` is not strictly necessary.
// TODO Handle terminating characters.
__attribute__ ((__nonnull__ (2)))
bool bmm_msg_preread(size_t* const ptr,
    struct bmm_msg_head const* const head) {
  // fprintf(stderr, "Prereading message %02x!\n", head->type);

  if (bmm_bit_test(head->flags, BMM_FBIT_BODY)) {
    if (bmm_bit_test(head->flags, BMM_FBIT_PREFIX)) {
      size_t const logsize = head->flags & 3;
      size_t const powsize = bmm_size_pow(2, logsize);

      // fprintf(stderr, "Log size %zu, power size %zu!\n", logsize, powsize);

      unsigned char buf[1 << 3];
      memset(buf, 0, sizeof buf); // Eh.
      switch (bmm_io_readin(buf, powsize)) {
        case BMM_IO_READ_ERROR:
        case BMM_IO_READ_EOF:
          BMM_ERR_FWARN(bmm_io_readin, "Failed to read size prefix");

          return false;
      }

      enum bmm_size_format const fmt =
        bmm_bit_test(head->flags, BMM_FBIT_INTLE) ? BMM_SIZE_FORMAT_LE :
        BMM_SIZE_FORMAT_BE;

      size_t bodysize;
      if (!bmm_size_from_buffer(&bodysize, buf,
            MAX(sizeof buf, sizeof bodysize), fmt)) {
        BMM_ERR_FWARN(bmm_size_from_buffer, "Message body too large");

        return false;
      }

      // fprintf(stderr, "Prefix has size %zu, contents '%02x %02x %02x %02x %02x %02x %02x %02x'.\n",
      //     bodysize, buf[0], buf[1], buf[2], buf[3], buf[4], buf[5], buf[6], buf[7]);

      if (ptr != NULL)
        *ptr = bodysize;
    } /* else {
      // TODO This would require an extra argument or more refactoring.
      unsigned char term;
      switch (bmm_io_readin(&term, 1)) {
        case BMM_IO_READ_ERROR:
        case BMM_IO_READ_EOF:
          BMM_ERR_FWARN(bmm_io_readin, "Failed to read terminator");

          return false;
      }

      // TODO This would be split into `bmm_msg_read` equivalent.
      // This is slow by design.
      for ever {
        unsigned char buf;
        switch (bmm_io_readin(&buf, 1)) {
          case BMM_IO_READ_ERROR:
          case BMM_IO_READ_EOF:
            BMM_ERR_FWARN(bmm_io_readin, "Failed to read body");

            return false;
        }

        if (buf == term)
          break;

        if (!bmm_io_writeout(&buf, 1)) {
          BMM_ERR_FWARN(bmm_io_writeout, "Failed to write body");

          return false;
        }
      }
    } */
  }

  return true;
}

__attribute__ ((__nonnull__ (1)))
bool bmm_msg_read(struct bmm_msg_head const* const head,
    void* const ptr, size_t const size) {
  size_t bodysize;
  if (!bmm_msg_preread(&bodysize, head))
    return false;

  if (!bmm_msg_precheck(size, bodysize))
    return false;

  switch (bmm_io_readin(ptr, bodysize)) {
    case BMM_IO_READ_ERROR:
    case BMM_IO_READ_EOF:
      BMM_ERR_FWARN(bmm_io_readin, "Failed to read body");

      return false;
  }

  return true;
}

// Assuming `head` has not been written...
__attribute__ ((__nonnull__))
bool bmm_msg_prewrite(struct bmm_msg_head const* const bad_head,
    size_t const size) {
  enum bmm_size_format const fmt =
    bmm_bit_test(bad_head->flags, BMM_FBIT_INTLE) ? BMM_SIZE_FORMAT_LE :
    BMM_SIZE_FORMAT_BE;

  // fprintf(stderr, "Prewriting message %02x, size %zu!\n", bad_head->type, size);

  size_t m;
  unsigned char buf[1 << 3];
  memset(buf, 0, sizeof buf); // Eh.
  if (!bmm_size_to_buffer(&m, buf, sizeof buf, size, fmt))
    return false;

  // fprintf(stderr, "Prefix has size %zu, contents '%02x %02x %02x %02x %02x %02x %02x %02x'.\n",
  //     m, buf[0], buf[1], buf[2], buf[3], buf[4], buf[5], buf[6], buf[7]);

  size_t const logsize = bmm_size_cilog(m, 2);
  size_t const powsize = bmm_size_pow(2, logsize);

  // fprintf(stderr, "Log size %zu, power size %zu!\n", logsize, powsize);

  struct bmm_msg_head head;
  head = *bad_head;
  head.flags = (head.flags & (~3)) | (logsize & 3);
  bmm_bit_pset(&head.flags, BMM_FBIT_BODY);
  bmm_bit_pset(&head.flags, BMM_FBIT_PREFIX);

  return bmm_io_writeout(&head, sizeof head) &&
    bmm_io_writeout(buf, powsize);
}

__attribute__ ((__nonnull__ (1)))
bool bmm_msg_write(struct bmm_msg_head const* const head,
    void const* const ptr, size_t const size) {
  if (!bmm_msg_prewrite(head, size))
    return false;

  return bmm_io_writeout(ptr, size);
}

bool bmm_msg_get(struct bmm_msg_head* const head, struct bmm_dem* const dem) {
  switch (bmm_io_readin(head, sizeof *head)) {
    case BMM_IO_READ_ERROR:
      BMM_ERR_FWARN(NULL, "Failed to read message header");

      return false;
    case BMM_IO_READ_EOF:
      return true;
  }

  struct bmm_dem_buf* const buf = bmm_dem_getbuf(dem);

  size_t bodysize;

  switch (head->type) {
    case BMM_MSG_NOP:
      if (!bmm_msg_read(head, NULL, 0))
        BMM_ERR_FWARN(NULL, "Failed to read nothing");
      break;
    case BMM_MSG_EKINE:
      // sizeof dem->istep + sizeof dem->est
      if (!bmm_msg_preread(&bodysize, head)) {
        BMM_ERR_FWARN(bmm_io_readin, "Failed to read estimators");

        return false;
      }
      switch (bmm_io_readin(&dem->istep, sizeof dem->istep)) {
        case BMM_IO_READ_ERROR:
        case BMM_IO_READ_EOF:
          BMM_ERR_FWARN(bmm_io_readin, "Failed to read estimators");

          return false;
      }
      switch (bmm_io_readin(&dem->est, sizeof dem->est)) {
        case BMM_IO_READ_ERROR:
        case BMM_IO_READ_EOF:
          BMM_ERR_FWARN(bmm_io_readin, "Failed to read estimators");

          return false;
      }
      break;
    case BMM_MSG_NPART:
      if (!bmm_msg_read(head, &buf->npart, sizeof buf->npart))
        BMM_ERR_FWARN(NULL, "Failed to read number of particles");
      break;
    case BMM_MSG_PARTS:
      // sizeof dem->istep + sizeof buf->parts + sizeof buf->partcs
      if (!bmm_msg_preread(&bodysize, head)) {
        BMM_ERR_FWARN(bmm_io_readin, "Failed to read particles");

        return false;
      }
      switch (bmm_io_readin(&dem->istep, sizeof dem->istep)) {
        case BMM_IO_READ_ERROR:
        case BMM_IO_READ_EOF:
          BMM_ERR_FWARN(bmm_io_readin, "Failed to read particles");

          return false;
      }
      switch (bmm_io_readin(&buf->parts, sizeof buf->parts)) {
        case BMM_IO_READ_ERROR:
        case BMM_IO_READ_EOF:
          BMM_ERR_FWARN(bmm_io_readin, "Failed to read particles");

          return false;
      }
      switch (bmm_io_readin(&buf->partcs, sizeof buf->partcs)) {
        case BMM_IO_READ_ERROR:
        case BMM_IO_READ_EOF:
          BMM_ERR_FWARN(bmm_io_readin, "Failed to read particles");

          return false;
      }
      break;
    case BMM_MSG_NEIGH:
      // sizeof buf->neigh + sizeof buf->links
      if (!bmm_msg_preread(&bodysize, head)) {
        BMM_ERR_FWARN(bmm_io_readin, "Failed to read neighbors");

        return false;
      }
      switch (bmm_io_readin(&buf->neigh, sizeof buf->neigh)) {
        case BMM_IO_READ_ERROR:
        case BMM_IO_READ_EOF:
          BMM_ERR_FWARN(bmm_io_readin, "Failed to read neighbors");

          return false;
      }
      switch (bmm_io_readin(&buf->links, sizeof buf->links)) {
        case BMM_IO_READ_ERROR:
        case BMM_IO_READ_EOF:
          BMM_ERR_FWARN(bmm_io_readin, "Failed to read neighbors");

          return false;
      }
      break;
    default:
      BMM_ERR_FWARN(NULL, "Unsupported message type");
  }

  return true;
}

void bmm_msg_put(struct bmm_msg_head const* const head,
    struct bmm_dem const* const dem) {
  struct bmm_dem_buf const* const buf = bmm_dem_getrbuf(dem);

  switch (head->type) {
    case BMM_MSG_NOP:
      if (!bmm_msg_write(head, NULL, 0))
        BMM_ERR_FWARN(NULL, "Failed to write nothing");
      break;
    case BMM_MSG_EKINE:
      if (!bmm_msg_prewrite(head, sizeof dem->istep + sizeof dem->est))
        BMM_ERR_FWARN(NULL, "Failed to write stuff");
      if (!bmm_io_writeout(&dem->istep, sizeof dem->istep))
        BMM_ERR_FWARN(NULL, "Failed to write particles");
      if (!bmm_io_writeout(&dem->est, sizeof dem->est))
        BMM_ERR_FWARN(NULL, "Failed to write estimators");
      break;
    case BMM_MSG_NPART:
      if (!bmm_msg_write(head, &buf->npart, sizeof buf->npart))
        BMM_ERR_FWARN(NULL, "Failed to write number of particles");
      break;
    case BMM_MSG_PARTS:
      if (!bmm_msg_prewrite(head,
            sizeof dem->istep + sizeof buf->parts + sizeof buf->partcs))
        BMM_ERR_FWARN(NULL, "Failed to write stuff");
      if (!bmm_io_writeout(&dem->istep, sizeof dem->istep))
        BMM_ERR_FWARN(NULL, "Failed to write particles");
      if (!bmm_io_writeout(&buf->parts, sizeof buf->parts))
        BMM_ERR_FWARN(NULL, "Failed to write particles");
      if (!bmm_io_writeout(&buf->partcs, sizeof buf->partcs))
        BMM_ERR_FWARN(NULL, "Failed to write particles");
      break;
    case BMM_MSG_NEIGH:
      if (!bmm_msg_prewrite(head, sizeof buf->neigh + sizeof buf->links))
        BMM_ERR_FWARN(NULL, "Failed to write stuff");
      if (!bmm_io_writeout(&buf->neigh, sizeof buf->neigh))
        BMM_ERR_FWARN(NULL, "Failed to write neighbors");
      if (!bmm_io_writeout(&buf->links, sizeof buf->links))
        BMM_ERR_FWARN(NULL, "Failed to write links");
      break;
    default:
      BMM_ERR_FWARN(NULL, "Unsupported message type");
  }

  if (bmm_bit_test(head->flags, BMM_FBIT_FLUSH))
    if (fflush(stdout) == EOF)
      BMM_ERR_FWARN(NULL, "Failed to flush output stream");
}
