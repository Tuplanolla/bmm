#include "bit.h"
#include "dem.h"
#include "err.h"
#include "io.h"
#include "msg.h"
#include <stdbool.h>
#include <stdio.h>

void bmm_defhead(struct bmm_msg_head* const head) {
  head->flags = 0;
  head->type = 0;
}

typedef size_t (* fproc)(void*, size_t);

static size_t slurp(void* const ptr, size_t const size) {
  return fread(ptr, size, 1, stdin);
}

static size_t ignore(__attribute__ ((__unused__)) void* const ptr,
    size_t const size) {
  return bmm_io_fastfw(stdin, size) / size;
}

static size_t barf(void const* const ptr, size_t const size) {
  return fwrite(ptr, size, 1, stdout);
}

bool bmm_msg_get(struct bmm_msg_head* const head, struct bmm_dem* const dem,
    bool (* const f)(struct bmm_msg_head const*, void*), void* const ptr) {
  if (slurp(head, sizeof *head) != 1) {
    if (feof(stdin) != 0)
      return false;

    BMM_ERR_FWARN(NULL, "Failed to read message header.");
  }

  bool const p = f(head, ptr);

  fproc proc = p ? slurp : ignore;

  switch (head->type) {
    case BMM_MSG_NOP:
      break;
    case BMM_MSG_NPART:
      if (proc(&dem->opts.npart, sizeof dem->opts.npart) != 1)
        BMM_ERR_FWARN(NULL, "Failed to read number of particles.");
      break;
    case BMM_MSG_PARTS:
      {
        struct bmm_dem_buf* const buf = bmm_dem_getwbuf(dem);
        if (proc(&buf->parts, sizeof buf->parts) != 1)
          BMM_ERR_FWARN(NULL, "Failed to read particles.");
      }
      break;
    default:
      BMM_ERR_FWARN(NULL, "Unsupported message type.");
  }

  return p;
}

void bmm_msg_put(struct bmm_msg_head const* const head,
    struct bmm_dem const* const dem) {
  if (barf(head, sizeof *head) != 1)
    BMM_ERR_FWARN(NULL, "Failed to write message header.");

  switch (head->type) {
    case BMM_MSG_NOP:
      break;
    case BMM_MSG_NPART:
      if (barf(&dem->opts.npart, sizeof dem->opts.npart) != 1)
        BMM_ERR_FWARN(NULL, "Failed to write number of particles.");
      break;
    case BMM_MSG_PARTS:
      {
        struct bmm_dem_buf const* const buf = bmm_dem_getrbuf(dem);
        if (barf(&buf->parts, sizeof buf->parts) != 1)
          BMM_ERR_FWARN(NULL, "Failed to write particles.");
      }
      break;
    default:
      BMM_ERR_FWARN(NULL, "Unsupported message type.");
  }

  if (bmm_bit_test(head->type, BMM_FBIT_FLUSH))
    if (fflush(stdout) == EOF)
      BMM_ERR_FWARN(NULL, "Failed to flush output stream.");
}
