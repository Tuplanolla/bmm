#include "bit.h"
#include "dem.h"
#include "err.h"
#include "filter.h"
#include "io.h"
#include "msg.h"
#include "size.h"
#include <stdio.h>

void bmm_filter_defopts(struct bmm_filter_opts* const opts) {
  for (size_t imsg = 0; imsg < BMM_MSG_MAX; ++imsg)
    opts->mask[imsg] = true;
}

void bmm_filter_def(struct bmm_filter* const filter,
    struct bmm_filter_opts const* const opts) {
  struct bmm_dem_opts defopts;
  bmm_dem_defopts(&defopts);
  bmm_dem_def(&filter->dem, &defopts);

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

    if (f(filter, &head)) {
      bmm_io_writeout(&head, sizeof head);

      if (bmm_bit_test(head.flags, BMM_FBIT_BODY)) {
        if (bmm_bit_test(head.flags, BMM_FBIT_PREFIX)) {
          size_t const presize = bmm_size_pow(2, head.flags & 3);

          unsigned char buf[1 << 3];
          switch (bmm_io_readin(buf, presize)) {
            case BMM_IO_READ_ERROR:
            case BMM_IO_READ_EOF:
              BMM_ERR_FWARN(bmm_io_readin, "Failed to read size prefix");

              return false;
          }

          enum bmm_size_format const fmt =
            bmm_bit_test(head.flags, BMM_FBIT_INTLE) ? BMM_SIZE_FORMAT_LE :
            BMM_SIZE_FORMAT_BE;

          size_t bodysize;
          if (!bmm_size_from_buffer(&bodysize, buf,
                MAX(sizeof buf, sizeof bodysize), fmt)) {
            BMM_ERR_FWARN(bmm_size_from_buffer, "Message body too large");

            return false;
          }

          if (!bmm_io_redirio(bodysize)) {
            BMM_ERR_FWARN(bmm_io_redirio, "Failed to redirect body");

            return false;
          }
        } else {
          unsigned char term;
          switch (bmm_io_readin(&term, 1)) {
            case BMM_IO_READ_ERROR:
            case BMM_IO_READ_EOF:
              BMM_ERR_FWARN(bmm_io_readin, "Failed to read terminator");

              return false;
          }

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
        }
      }

      if (bmm_bit_test(head.flags, BMM_FBIT_FLUSH))
        if (fflush(stdout) == EOF)
          return false;
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
