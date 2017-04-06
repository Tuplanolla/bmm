#include <signal.h>
#include <stdio.h>

#include "bit.h"
#include "dem.h"
#include "err.h"
#include "filter.h"
#include "io.h"
#include "msg.h"
#include "sig.h"
#include "size.h"

void bmm_filter_opts_def(struct bmm_filter_opts* const opts) {
  for (size_t imsg = 0; imsg < BMM_MSG_MAX; ++imsg)
    opts->mask[imsg] = false;
}

void bmm_filter_def(struct bmm_filter* const filter,
    struct bmm_filter_opts const* const opts) {
  filter->opts = *opts;
  filter->passed = 0;
  filter->stopped = 0;
}

static bool f(struct bmm_filter const* const filter,
    struct bmm_msg_head const* const head) {
  return filter->opts.mask[head->type];
}

bool bmm_filter_run(struct bmm_filter* const filter) {
  int const sigs[] = {SIGINT, SIGQUIT, SIGTERM, SIGPIPE};
  if (bmm_sig_register(sigs, sizeof sigs / sizeof *sigs) != SIZE_MAX) {
    BMM_ERR_WARN(bmm_sig_register);

    return false;
  }

  for ever {
    int signum;
    if (bmm_sig_use(&signum))
      switch (signum) {
        case SIGINT:
        case SIGQUIT:
        case SIGTERM:
        case SIGPIPE:
          BMM_ERR_FWARN(NULL, "Filtering interrupted");

          return false;
      }

    struct bmm_msg_head head;
    switch (bmm_io_readin(&head, sizeof head)) {
      case BMM_IO_READ_ERROR:
        BMM_ERR_FWARN(bmm_io_readin, "Failed to read header");

        return false;
      case BMM_IO_READ_EOF:
        return true;
    }

    size_t bodysize;
    if (!bmm_msg_preread(&bodysize, &head)) {
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

      ++filter->passed;
    } else {
      if (!bmm_io_fastfwin(bodysize)) {
        BMM_ERR_FWARN(bmm_io_fastfwin, "Failed to fast-forward body");

        return false;
      }

      ++filter->stopped;
    }
  }

  return true;
}

bool bmm_filter_run_with(struct bmm_filter_opts const* const opts) {
  struct bmm_filter filter;
  bmm_filter_def(&filter, opts);

  return bmm_filter_run(&filter);
}
