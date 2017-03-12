#include "bits.h"
#include "dem.h"
#include "errors.h"
#include "msgs.h"
#include <stdio.h>

__attribute__ ((__nonnull__))
void bmm_defhead(struct bmm_head* const head) {
  head->flags = 0;
  head->type = 0;
}

// TODO Check sizes etc.
void bmm_getmsg(FILE* const stream, struct bmm_head* const head,
    struct bmm_state* const state) {
  if (fread(head, sizeof *head, 1, stream) != 1)
    bmm_error("Failed to read message header.");

  switch (head->type) {
    case BMM_MSG_NDIM:
      if (fread(&state->opts.ndim, sizeof state->opts.ndim, 1, stream) != 1)
        bmm_error("Failed to read number of dimensions.");
      break;
    case BMM_MSG_PARTS:
      if (fread(&state->parts, sizeof state->parts, 1, stream) != 1)
        bmm_error("Failed to read particles.");
      break;
    default:
      bmm_error("Unsupported message type.");
  }
}

void bmm_putmsg(FILE* const stream, struct bmm_head const* const head,
    struct bmm_state const* const state) {
  if (fwrite(head, sizeof *head, 1, stream) != 1)
    bmm_error("Failed to write message header.");

  switch (head->type) {
    case BMM_MSG_NDIM:
      if (fwrite(&state->opts.ndim, sizeof state->opts.ndim, 1, stream) != 1)
        bmm_error("Failed to write number of dimensions.");
      break;
    case BMM_MSG_PARTS:
      if (fwrite(&state->parts, sizeof state->parts, 1, stream) != 1)
        bmm_error("Failed to write particles.");
      break;
    default:
      bmm_error("Unsupported message type.");
  }

  if (bmm_testb(head->type, BMM_FBIT_FLUSH))
    if (fflush(stream) == EOF)
      bmm_error("Failed to flush output stream");
}
