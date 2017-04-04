#include "bit.h"
#include "dem.h"
#include "err.h"
#include "filter.h"
#include "io.h"
#include "msg.h"
#include "size.h"
#include <stdint.h> // TODO Remove temporary SIZE_MAX.
#include <stdio.h>

void bmm_filter_defopts(struct bmm_filter_opts* const opts) {
  for (size_t imsg = 0; imsg < BMM_MSG_MAX; ++imsg)
    opts->mask[imsg] = true;
}

void bmm_filter_def(struct bmm_filter* const filter,
    struct bmm_filter_opts const* const opts) {
  struct bmm_dem_opts defopts;
  bmm_dem_defopts(&defopts);

  filter->opts = *opts;
}

static bool f(struct bmm_filter const* const filter,
    struct bmm_msg_head const* const head) {
  return filter->opts.mask[head->type];
}

bool bmm_filter_run_with(struct bmm_filter* const filter) {
  for ever {
    struct bmm_msg_head head;
    switch (bmm_io_readin(&head, sizeof head)) {
      case BMM_IO_READ_ERROR:
        BMM_ERR_FWARN(bmm_io_readin, "Failed to read header");

        return false;
      case BMM_IO_READ_EOF:
        return true;
    }

    size_t bodysize;
    if (!bmm_msg_preread(&bodysize, &head, SIZE_MAX)) {
      BMM_ERR_FWARN(bmm_io_readin, "Failed to read prefix (and do stuff)");

      return false;
    }

    if (f(filter, &head)) {
      if (!bmm_msg_prewrite(&head, bodysize))
        BMM_ERR_FWARN(NULL, "Failed to write stuff");

      if (!bmm_io_redirio(bodysize)) {
        BMM_ERR_FWARN(bmm_io_redirio, "Failed to redirect body");

        return false;
      }

      if (bmm_bit_test(head.flags, BMM_FBIT_FLUSH))
        if (fflush(stdout) == EOF)
          return false;
    } else {
      if (!bmm_io_fastfwin(bodysize)) {
        BMM_ERR_FWARN(bmm_io_fastfwin, "Failed to fast-forward body");

        return false;
      }
    }
  }

  return true;
}

bool bmm_filter_run(struct bmm_filter_opts const* const opts) {
  struct bmm_filter* const filter = malloc(sizeof *filter);
  if (filter == NULL) {
    BMM_ERR_WARN(malloc);

    return false;
  }

  bmm_filter_def(filter, opts);
  bool const result = bmm_filter_run_with(filter);

  free(filter);

  return result;
}
