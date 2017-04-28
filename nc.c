#include <netcdf.h>
#include <signal.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "ext.h"
#include "io.h"
#include "msg.h"
#include "nc.h"
#include "sig.h"
#include "tle.h"

void bmm_nc_opts_def(struct bmm_nc_opts* const opts) {
  opts->conv = BMM_NC_CONV_AMBER;

  char const path[] = "bmm.nc";
  static_assert(sizeof opts->path >= sizeof path, "Path too long");
  (void) memcpy(opts->path, path, sizeof path);

  opts->i = false;
  opts->r = true;
  opts->q = true;
  opts->v = false;
  opts->f = false;
}

void bmm_nc_def(struct bmm_nc* const nc,
    struct bmm_nc_opts const* const opts) {
  nc->opts = *opts;
  nc->stream = NULL;
  nc->npart = 0;
}

static enum bmm_io_read msg_read(void* buf, size_t const n,
    __attribute__ ((__unused__)) void* const ptr) {
  return bmm_io_readin(buf, n);
}

enum bmm_io_read bmm_nc_step(struct bmm_nc* const nc) {
  struct bmm_msg_spec spec;
  switch (bmm_msg_spec_read(&spec, msg_read, NULL)) {
    case BMM_IO_READ_ERROR:
      BMM_TLE_EXTS(BMM_TLE_NUM_IO, "Failed to read message header");

      return BMM_IO_READ_ERROR;
    case BMM_IO_READ_EOF:
      return BMM_IO_READ_EOF;
  }

  if (spec.endy != bmm_endy_get()) {
    BMM_TLE_EXTS(BMM_TLE_NUM_UNIMPL, "Unsupported endianness");

    return BMM_IO_READ_ERROR;
  }

  if (spec.tag != BMM_MSG_TAG_SP) {
    BMM_TLE_EXTS(BMM_TLE_NUM_UNIMPL, "Unsupported tag");

    return BMM_IO_READ_ERROR;
  }

  enum bmm_msg_num num;
  switch (bmm_msg_num_read(&num, msg_read, NULL)) {
    case BMM_IO_READ_ERROR:
    case BMM_IO_READ_EOF:
      BMM_TLE_EXTS(BMM_TLE_NUM_IO, "Failed to read message number");

      return BMM_IO_READ_ERROR;
  }

  switch (num) {
    case BMM_MSG_NUM_NPART:
      ; // TODO Actually do some work.
  }

  return BMM_IO_READ_SUCCESS;
}

bool bmm_nc_run(struct bmm_nc* const nc) {
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
          BMM_TLE_EXTS(BMM_TLE_NUM_ASYNC, "Exporting interrupted");

          return false;
      }

    switch (bmm_nc_step(nc)) {
      case BMM_IO_READ_ERROR:
        return false;
      case BMM_IO_READ_EOF:
        return true;
    }
  }

  return true;
}

bool bmm_nc_run_with(struct bmm_nc_opts const* const opts) {
  struct bmm_nc nc;
  bmm_nc_def(&nc, opts);

  return bmm_nc_run(&nc);
}
