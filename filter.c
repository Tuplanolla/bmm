#include <signal.h>
#include <stdio.h>

#include "bit.h"
#include "filter.h"
#include "io.h"
#include "msg.h"
#include "sig.h"
#include "size.h"
#include "tle.h"

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

static bool pass(struct bmm_filter const* const filter,
    enum bmm_msg_type const type) {
  return filter->opts.mask[(size_t) type];
}

static enum bmm_io_read msg_read(void* buf, size_t const n,
    __attribute__ ((__unused__)) void* const ptr) {
  return bmm_io_readin(buf, n);
}

static bool msg_write(void const* buf, size_t const n,
    __attribute__ ((__unused__)) void* const ptr) {
  return bmm_io_writeout(buf, n);
}

bool bmm_filter_run(struct bmm_filter* const filter) {
  int const sigs[] = {SIGINT, SIGQUIT, SIGTERM, SIGPIPE};
  if (bmm_sig_register(sigs, sizeof sigs / sizeof *sigs) != SIZE_MAX) {
    BMM_TLE_STDS();

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
          BMM_TLE_EXTS(BMM_TLE_ASYNC, "Filtering interrupted");

          return false;
      }

    struct bmm_msg_spec spec;
    switch (bmm_msg_spec_read(&spec, msg_read, NULL)) {
      case BMM_IO_READ_ERROR:
        BMM_TLE_EXTS(BMM_TLE_IO, "Failed to read header");

        return false;
      case BMM_IO_READ_EOF:
        return true;
    }

    enum bmm_msg_type type;
    switch (bmm_msg_type_read(&type, msg_read, NULL)) {
      case BMM_IO_READ_ERROR:
      case BMM_IO_READ_EOF:
        BMM_TLE_EXTS(BMM_TLE_IO, "Failed to read type");

        return false;
    }

    if (pass(filter, type)) {
      if (!bmm_msg_spec_write(&spec, msg_write, NULL) ||
          !bmm_msg_type_write(&type, msg_write, NULL)) {
        BMM_TLE_EXTS(BMM_TLE_IO, "Failed to write stuff");

        return false;
      }

      switch (spec.tag) {
        case BMM_MSG_TAG_SP:
          if (!bmm_io_redirio(spec.msg.size - BMM_MSG_TYPESIZE)) {
            BMM_TLE_EXTS(BMM_TLE_IO, "Failed to redirect body");

            return false;
          }

          break;
        case BMM_MSG_TAG_LT:
          BMM_TLE_EXTS(BMM_TLE_UNIMPL, "Not implemented (yet?)");

          return false;
      }

      ++filter->passed;
    } else {
      switch (spec.tag) {
        case BMM_MSG_TAG_SP:
          if (!bmm_io_fastfwin(spec.msg.size - BMM_MSG_TYPESIZE)) {
            BMM_TLE_EXTS(BMM_TLE_IO, "Failed to fast-forward body");

            return false;
          }

          break;
        case BMM_MSG_TAG_LT:
          BMM_TLE_EXTS(BMM_TLE_UNIMPL, "Not implemented (yet?)");

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
