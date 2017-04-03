#include "bit.h"
#include "dem.h"
#include "err.h"
#include "io.h"
#include "msg.h"
#include "size.h"
#include <stdbool.h>
#include <stdio.h>

void bmm_defhead(struct bmm_msg_head* const head) {
  head->flags = 0;
  head->type = 0;
}

// Assuming `head` has already been read...
__attribute__ ((__nonnull__))
bool bmm_msg_read(struct bmm_msg_head const* const head,
    void* const ptr, size_t const size) {
  if (bmm_bit_test(head->flags, BMM_FBIT_BODY)) {
    if (bmm_bit_test(head->flags, BMM_FBIT_PREFIX)) {
      size_t const presize = bmm_size_pow(2, head->flags & 3);

      unsigned char buf[1 << 3];
      switch (bmm_io_readin(buf, presize)) {
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

      if (size < bodysize) {
        BMM_ERR_FWARN(bmm_size_from_buffer, "Message body would overflow");

        return false;
      }

      if (size != bodysize)
        BMM_ERR_FWARN(bmm_size_from_buffer, "Message body size mismatch");

    // TODO Factor out? This would allow receiving compound items separately.
      switch (bmm_io_readin(ptr, bodysize)) {
        case BMM_IO_READ_ERROR:
        case BMM_IO_READ_EOF:
          BMM_ERR_FWARN(bmm_io_readin, "Failed to read body");

          return false;
      }
    }
  }

  return true;
}

// Assuming `head` has already been written...
__attribute__ ((__nonnull__))
bool bmm_msg_write(struct bmm_msg_head const* const head,
    void const* const ptr, size_t const size) {
  enum bmm_size_format const fmt =
    bmm_bit_test(head->flags, BMM_FBIT_INTLE) ? BMM_SIZE_FORMAT_LE :
    BMM_SIZE_FORMAT_BE;

  size_t m;
  unsigned char buf[1 << 3];
  if (!bmm_size_to_buffer(&m, buf, sizeof buf, size, fmt))
    return false;

  return bmm_io_writeout(buf, m) &&
    // TODO Factor out? This would allow sending disjoint items together.
    bmm_io_writeout(ptr, size);
}

bool bmm_msg_get(struct bmm_msg_head* const head, struct bmm_dem* const dem,
    bool (* const f)(struct bmm_msg_head const*, void*), void* const ptr) {
  if (BMM_IO_READ_ERROR == bmm_io_readin(head, sizeof *head)) {
    if (feof(stdin) != 0)
      return false;

    BMM_ERR_FWARN(NULL, "Failed to read message header");
  }

  bool const p = f(head, ptr);

  enum bmm_io_read (* proc)(void*, size_t) = p ? bmm_io_readin : bmm_io_fastfwin;

  struct bmm_dem_buf* const buf = bmm_dem_getbuf(dem);

  switch (head->type) {
    case BMM_MSG_NOP:
      break;
    case BMM_MSG_EKINE:
      if (BMM_IO_READ_ERROR == proc(&dem->istep, sizeof dem->istep))
        BMM_ERR_FWARN(NULL, "Failed to write particles");
      if (BMM_IO_READ_ERROR == proc(&dem->est, sizeof dem->est))
        BMM_ERR_FWARN(NULL, "Failed to read estimators");
      break;
    case BMM_MSG_NPART:
      if (BMM_IO_READ_ERROR == proc(&buf->npart, sizeof buf->npart))
        BMM_ERR_FWARN(NULL, "Failed to read number of particles");
      break;
    case BMM_MSG_PARTS:
      if (BMM_IO_READ_ERROR == proc(&dem->istep, sizeof dem->istep))
        BMM_ERR_FWARN(NULL, "Failed to read particles");
      if (BMM_IO_READ_ERROR == proc(&buf->parts, sizeof buf->parts))
        BMM_ERR_FWARN(NULL, "Failed to read particles");
      if (BMM_IO_READ_ERROR == proc(&buf->partcs, sizeof buf->partcs))
        BMM_ERR_FWARN(NULL, "Failed to read particles");
      break;
    case BMM_MSG_NEIGH:
      if (BMM_IO_READ_ERROR == proc(&buf->neigh, sizeof buf->neigh))
        BMM_ERR_FWARN(NULL, "Failed to read neighbors");
      if (BMM_IO_READ_ERROR == proc(&buf->links, sizeof buf->links))
        BMM_ERR_FWARN(NULL, "Failed to read links");
      break;
    default:
      BMM_ERR_FWARN(NULL, "Unsupported message type");
  }

  return p;
}

void bmm_msg_put(struct bmm_msg_head const* const head,
    struct bmm_dem const* const dem) {
  if (!bmm_io_writeout(head, sizeof *head))
    BMM_ERR_FWARN(NULL, "Failed to write message header");

  struct bmm_dem_buf const* const buf = bmm_dem_getrbuf(dem);

  switch (head->type) {
    case BMM_MSG_NOP:
      break;
    case BMM_MSG_EKINE:
      if (!bmm_io_writeout(&dem->istep, sizeof dem->istep))
        BMM_ERR_FWARN(NULL, "Failed to write particles");
      if (!bmm_io_writeout(&dem->est, sizeof dem->est))
        BMM_ERR_FWARN(NULL, "Failed to write estimators");
      break;
    case BMM_MSG_NPART:
      if (!bmm_io_writeout(&buf->npart, sizeof buf->npart))
        BMM_ERR_FWARN(NULL, "Failed to write number of particles");
      break;
    case BMM_MSG_PARTS:
      if (!bmm_io_writeout(&dem->istep, sizeof dem->istep))
        BMM_ERR_FWARN(NULL, "Failed to write particles");
      if (!bmm_io_writeout(&buf->parts, sizeof buf->parts))
        BMM_ERR_FWARN(NULL, "Failed to write particles");
      if (!bmm_io_writeout(&buf->partcs, sizeof buf->partcs))
        BMM_ERR_FWARN(NULL, "Failed to write particles");
      break;
    case BMM_MSG_NEIGH:
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
