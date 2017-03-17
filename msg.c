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

typedef size_t (* fproc)(void*, size_t, size_t, FILE*);

static size_t fignore(void* const ptr,
    size_t const size, size_t const nmemb, FILE* const stream) {
  return bmm_io_fastfw(size, nmemb, stream);
}

// TODO Check sizes etc.
bool bmm_msg_get(struct bmm_msg_head* const head, struct bmm_dem* const dem,
    bool (* const f)(struct bmm_msg_head const*, void*), void* const ptr) {
  if (fread(head, sizeof *head, 1, stdin) != 1) {
    if (feof(stdin) != 0)
      return false;

    BMM_ERR_WARN("Failed to read message header.");
  }

  bool const p = f(head, ptr);

  fproc fproc = p ? fread : fignore;

  switch (head->type) {
    case BMM_MSG_NOP:
      break;
    case BMM_MSG_NDIM:
      if (fproc(&dem->opts.ndim, sizeof dem->opts.ndim, 1, stdin) != 1)
        BMM_ERR_WARN("Failed to read number of dimensions.");
      break;
    case BMM_MSG_PARTS:
      if (fproc(&dem->parts, sizeof dem->parts, 1, stdin) != 1)
        BMM_ERR_WARN("Failed to read particles.");
      break;
    default:
      BMM_ERR_WARN("Unsupported message type.");
  }

  return p;
}

void bmm_msg_put(struct bmm_msg_head const* const head,
    struct bmm_dem const* const dem) {
  if (fwrite(head, sizeof *head, 1, stdout) != 1)
    BMM_ERR_WARN("Failed to write message header.");

  switch (head->type) {
    case BMM_MSG_NOP:
      break;
    case BMM_MSG_NDIM:
      if (fwrite(&dem->opts.ndim, sizeof dem->opts.ndim, 1, stdout) != 1)
        BMM_ERR_WARN("Failed to write number of dimensions.");
      break;
    case BMM_MSG_PARTS:
      if (fwrite(&dem->parts, sizeof dem->parts, 1, stdout) != 1)
        BMM_ERR_WARN("Failed to write particles.");
      break;
    default:
      BMM_ERR_WARN("Unsupported message type.");
  }

  if (bmm_bit_test(head->type, BMM_FBIT_FLUSH))
    if (fflush(stdout) == EOF)
      BMM_ERR_WARN("Failed to flush output stream.");
}
