#include "bit.h"
#include "dem.h"
#include "err.h"
#include "msg.h"
#include <stdbool.h>
#include <stdio.h>

void bmm_defhead(struct bmm_head* const head) {
  head->flags = 0;
  head->type = 0;
}

// TODO Check sizes etc.
bool bmm_msg_get(struct bmm_head* const head,
    struct bmm_dem* const dem) {
  if (fread(head, sizeof *head, 1, stdin) != 1) {
    if (feof(stdin) != 0)
      return false;

    BMM_ERR_WARN("Failed to read message header.");
  }

  switch (head->type) {
    case BMM_MSG_NDIM:
      if (fread(&dem->opts.ndim, sizeof dem->opts.ndim, 1, stdin) != 1)
        BMM_ERR_WARN("Failed to read number of dimensions.");
      break;
    case BMM_MSG_PARTS:
      if (fread(&dem->parts, sizeof dem->parts, 1, stdin) != 1)
        BMM_ERR_WARN("Failed to read particles.");
      break;
    default:
      BMM_ERR_WARN("Unsupported message type.");
  }

  return true;
}

void bmm_msg_put(struct bmm_head const* const head,
    struct bmm_dem const* const dem) {
  if (fwrite(head, sizeof *head, 1, stdout) != 1)
    BMM_ERR_WARN("Failed to write message header.");

  switch (head->type) {
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
