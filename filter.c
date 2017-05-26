#include <signal.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>

#include "conf.h"
#include "filter.h"
#include "fp.h"
#include "io.h"
#include "msg.h"
#include "sig.h"
#include "tle.h"

void bmm_filter_opts_def(struct bmm_filter_opts* const opts) {
  opts->verbose = false;

  for (size_t imsg = 0; imsg < BMM_MMSG; ++imsg)
    opts->mask[imsg] = false;
}

void bmm_filter_def(struct bmm_filter* const filter,
    struct bmm_filter_opts const* const opts) {
  filter->opts = *opts;
  filter->passed = 0;
  filter->stopped = 0;
}

static bool pass(struct bmm_filter const* const filter,
    enum bmm_msg_num const num) {
  return filter->opts.mask[(size_t) num];
}

static enum bmm_io_read msg_read(void* buf, size_t const n,
    __attribute__ ((__unused__)) void* const ptr) {
  return bmm_io_readin(buf, n);
}

static bool msg_write(void const* buf, size_t const n,
    __attribute__ ((__unused__)) void* const ptr) {
  return bmm_io_writeout(buf, n);
}

enum bmm_io_read bmm_filter_step(struct bmm_filter* const filter) {
  struct bmm_msg_spec spec;
  switch (bmm_msg_spec_read(&spec, msg_read, NULL)) {
    case BMM_IO_READ_ERROR:
      return BMM_IO_READ_ERROR;
    case BMM_IO_READ_EOF:
      return BMM_IO_READ_EOF;
  }

  enum bmm_msg_num num;
  switch (bmm_msg_num_read(&num, msg_read, NULL)) {
    case BMM_IO_READ_EOF:
      BMM_TLE_EXTS(BMM_TLE_NUM_IO, "Unexpected end");
    case BMM_IO_READ_ERROR:
      return BMM_IO_READ_ERROR;
  }

  if (pass(filter, num)) {
    if (!bmm_msg_spec_write(&spec, msg_write, NULL) ||
        !bmm_msg_num_write(&num, msg_write, NULL))
      return BMM_IO_READ_ERROR;

    switch (spec.tag) {
      case BMM_MSG_TAG_SP:
        if (!bmm_io_redirio(spec.msg.size - BMM_MSG_NUMSIZE))
          return BMM_IO_READ_ERROR;

        break;
      case BMM_MSG_TAG_LT:
        BMM_TLE_EXTS(BMM_TLE_NUM_UNIMPL, "Not implemented");

        return BMM_IO_READ_ERROR;
    }

    ++filter->passed;
  } else {
    switch (spec.tag) {
      case BMM_MSG_TAG_SP:
        switch (bmm_io_fastfwin(spec.msg.size - BMM_MSG_NUMSIZE)) {
          case BMM_IO_READ_EOF:
            BMM_TLE_EXTS(BMM_TLE_NUM_IO, "Unexpected end");
          case BMM_IO_READ_ERROR:
            return BMM_IO_READ_ERROR;
        }

        break;
      case BMM_MSG_TAG_LT:
        BMM_TLE_EXTS(BMM_TLE_NUM_UNIMPL, "Not implemented");

        return BMM_IO_READ_ERROR;
    }

    ++filter->stopped;
  }

  return BMM_IO_READ_SUCCESS;
}

bool bmm_filter_report(struct bmm_filter const* const filter) {
  if (filter->opts.verbose) {
    double const passed = filter->passed;
    double const stopped = filter->stopped;
    double const total = passed + stopped;

    if (fprintf(stderr, "Passed: %zu\n" "Stopped: %zu\n" "Ratio: %.f %%\n",
          filter->passed, filter->stopped, bmm_fp_percent(passed, total)) < 0)
      return false;
  }

  return true;
}

static bool bmm_filter_run_(struct bmm_filter* const filter) {
  int const sigs[] = {SIGINT, SIGQUIT, SIGTERM, SIGPIPE};
  if (bmm_sig_register(sigs, nmembof(sigs)) != SIZE_MAX) {
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
          BMM_TLE_EXTS(BMM_TLE_NUM_ASYNC, "Interrupted");

          return false;
      }

    switch (bmm_filter_step(filter)) {
      case BMM_IO_READ_ERROR:
        return false;
      case BMM_IO_READ_EOF:
        return true;
    }
  }

  return true;
}

bool bmm_filter_run(struct bmm_filter* const filter) {
  bool const run = bmm_filter_run_(filter);
  bool const report = bmm_filter_report(filter);

  return run && report;
}

bool bmm_filter_run_with(struct bmm_filter_opts const* const opts) {
  struct bmm_filter filter;
  bmm_filter_def(&filter, opts);

  return bmm_filter_run(&filter);
}
