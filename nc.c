#include <netcdf.h>
#include <signal.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "conf.h"
// TODO Undepend.
#include "dem.h"
#include "ext.h"
#include "io.h"
#include "math.h"
#include "msg.h"
#include "nc.h"
#include "sig.h"
#include "tle.h"

void bmm_nc_opts_def(struct bmm_nc_opts* const opts) {
  opts->conv = BMM_NC_CONV_AMBER;
  opts->path = "bmm.nc";
  opts->i = false;
  opts->r = true;
  opts->q = true;
  opts->v = false;
  opts->f = false;
}

// Must have `NDIM == 3` for OVITO.
#define NDIM 3

void bmm_nc_def(struct bmm_nc* const nc,
    struct bmm_nc_opts const* const opts) {
  nc->opts = *opts;
  // TODO This.
  // nc->npart = 0;
  nc->npart = 1;
  nc->iframe = 0;
}

static enum bmm_io_read msg_read(void* buf, size_t const n,
    __attribute__ ((__unused__)) void* const ptr) {
  return bmm_io_readin(buf, n);
}

bool bmm_nc_close(struct bmm_nc* const nc) {
  int nerr;

  nerr = nc_close(nc->ncid);
  if (nerr != NC_NOERR) {
    BMM_TLE_EXTS(BMM_TLE_NUM_NC, "%s", nc_strerror(nerr));

    return false;
  }

  return true;
}

bool bmm_nc_open(struct bmm_nc* const nc) {
  int nerr;

  nerr = nc_create(nc->opts.path, NC_CLOBBER | NC_64BIT_OFFSET, &nc->ncid);
  if (nerr != NC_NOERR) {
    BMM_TLE_EXTS(BMM_TLE_NUM_NC, "%s", nc_strerror(nerr));

    return false;
  }

  char const convs[] = "AMBER";
  nerr = nc_put_att_text(nc->ncid, NC_GLOBAL,
      "Conventions", strlen(convs), convs);
  if (nerr != NC_NOERR) {
    BMM_TLE_EXTS(BMM_TLE_NUM_NC, "%s", nc_strerror(nerr));

    return false;
  }

  char const convver[] = "1.0";
  nerr = nc_put_att_text(nc->ncid, NC_GLOBAL,
      "ConventionVersion", strlen(convver), convver);
  if (nerr != NC_NOERR) {
    BMM_TLE_EXTS(BMM_TLE_NUM_NC, "%s", nc_strerror(nerr));

    return false;
  }

  char const prog[] = BMM_PROGID;
  nerr = nc_put_att_text(nc->ncid, NC_GLOBAL,
      "program", strlen(prog), prog);
  if (nerr != NC_NOERR) {
    BMM_TLE_EXTS(BMM_TLE_NUM_NC, "%s", nc_strerror(nerr));

    return false;
  }

  char const progver[] = BMM_VERSION;
  nerr = nc_put_att_text(nc->ncid, NC_GLOBAL,
      "programVersion", strlen(progver), progver);
  if (nerr != NC_NOERR) {
    BMM_TLE_EXTS(BMM_TLE_NUM_NC, "%s", nc_strerror(nerr));

    return false;
  }

  nerr = nc_def_dim(nc->ncid, "frame", NC_UNLIMITED, &nc->id_frame);
  if (nerr != NC_NOERR) {
    BMM_TLE_EXTS(BMM_TLE_NUM_NC, "%s", nc_strerror(nerr));

    return false;
  }

  nerr = nc_def_dim(nc->ncid, "spatial", NDIM, &nc->id_spatial);
  if (nerr != NC_NOERR) {
    BMM_TLE_EXTS(BMM_TLE_NUM_NC, "%s", nc_strerror(nerr));

    return false;
  }

  nerr = nc_def_dim(nc->ncid, "atom", BMM_NPART, &nc->id_atom);
  if (nerr != NC_NOERR) {
    BMM_TLE_EXTS(BMM_TLE_NUM_NC, "%s", nc_strerror(nerr));

    return false;
  }

  nerr = nc_def_dim(nc->ncid, "cell_spatial", NDIM, &nc->id_cspatial);
  if (nerr != NC_NOERR) {
    BMM_TLE_EXTS(BMM_TLE_NUM_NC, "%s", nc_strerror(nerr));

    return false;
  }

  nerr = nc_def_dim(nc->ncid, "label", NDIM, &nc->id_label);
  if (nerr != NC_NOERR) {
    BMM_TLE_EXTS(BMM_TLE_NUM_NC, "%s", nc_strerror(nerr));

    return false;
  }

  int dimids[3];

  dimids[0] = nc->id_spatial;
  nerr = nc_def_var(nc->ncid, "spatial", NC_CHAR, 1, dimids, &nc->varid_spatial);
  if (nerr != NC_NOERR) {
    BMM_TLE_EXTS(BMM_TLE_NUM_NC, "%s", nc_strerror(nerr));

    return false;
  }

  dimids[0] = nc->id_cspatial;
  nerr = nc_def_var(nc->ncid, "cell_spatial", NC_CHAR, 1, dimids, &nc->varid_cspatial);
  if (nerr != NC_NOERR) {
    BMM_TLE_EXTS(BMM_TLE_NUM_NC, "%s", nc_strerror(nerr));

    return false;
  }

  dimids[0] = nc->id_frame;
  nerr = nc_def_var(nc->ncid, "time", NC_FLOAT, 1, dimids, &nc->varid_time);
  if (nerr != NC_NOERR) {
    BMM_TLE_EXTS(BMM_TLE_NUM_NC, "%s", nc_strerror(nerr));

    return false;
  }

  char const time[] = "second";
  nerr = nc_put_att_text(nc->ncid, nc->varid_time,
      "units", strlen(time), time);
  if (nerr != NC_NOERR) {
    BMM_TLE_EXTS(BMM_TLE_NUM_NC, "%s", nc_strerror(nerr));

    return false;
  }

  dimids[0] = nc->id_frame;
  dimids[1] = nc->id_atom;
  dimids[2] = nc->id_spatial;
  nerr = nc_def_var(nc->ncid,
      "coordinates", NC_FLOAT, 3, dimids, &nc->varid_coords);
  if (nerr != NC_NOERR) {
    BMM_TLE_EXTS(BMM_TLE_NUM_NC, "%s", nc_strerror(nerr));

    return false;
  }

  char const coords[] = "meter";
  nerr = nc_put_att_text(nc->ncid, nc->varid_coords,
      "units", strlen(coords), coords);
  if (nerr != NC_NOERR) {
    BMM_TLE_EXTS(BMM_TLE_NUM_NC, "%s", nc_strerror(nerr));

    return false;
  }

  dimids[0] = nc->id_frame;
  dimids[1] = nc->id_cspatial;
  nerr = nc_def_var(nc->ncid,
      "cell_lengths", NC_FLOAT, 2, dimids, &nc->varid_clens);
  if (nerr != NC_NOERR) {
    BMM_TLE_EXTS(BMM_TLE_NUM_NC, "%s", nc_strerror(nerr));

    return false;
  }

  char const clens[] = "meter";
  nerr = nc_put_att_text(nc->ncid, nc->varid_clens,
      "units", strlen(clens), clens);
  if (nerr != NC_NOERR) {
    BMM_TLE_EXTS(BMM_TLE_NUM_NC, "%s", nc_strerror(nerr));

    return false;
  }

  dimids[0] = nc->id_frame;
  dimids[1] = nc->id_atom;
  // Must have `name = "radius"` for OVITO.
  nerr = nc_def_var(nc->ncid,
      "radius", NC_FLOAT, 2, dimids, &nc->varid_radii);
  if (nerr != NC_NOERR) {
    BMM_TLE_EXTS(BMM_TLE_NUM_NC, "%s", nc_strerror(nerr));

    return false;
  }

  char const radii[] = "meter";
  nerr = nc_put_att_text(nc->ncid, nc->varid_radii,
      "units", strlen(radii), radii);
  if (nerr != NC_NOERR) {
    BMM_TLE_EXTS(BMM_TLE_NUM_NC, "%s", nc_strerror(nerr));

    return false;
  }

  nerr = nc_enddef(nc->ncid);
  if (nerr != NC_NOERR) {
    BMM_TLE_EXTS(BMM_TLE_NUM_NC, "%s", nc_strerror(nerr));

    return false;
  }

  nerr = nc_put_var_text(nc->ncid, nc->varid_spatial, "xyz");
  if (nerr != NC_NOERR) {
    BMM_TLE_EXTS(BMM_TLE_NUM_NC, "%s", nc_strerror(nerr));

    return false;
  }

  nerr = nc_put_var_text(nc->ncid, nc->varid_cspatial, "abc");
  if (nerr != NC_NOERR) {
    BMM_TLE_EXTS(BMM_TLE_NUM_NC, "%s", nc_strerror(nerr));

    return false;
  }

  return true;
}

enum bmm_io_read bmm_nc_step(struct bmm_nc* const nc) {
  struct bmm_msg_spec spec;
  switch (bmm_msg_spec_read(&spec, msg_read, NULL)) {
    case BMM_IO_READ_ERROR:
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
    case BMM_IO_READ_EOF:
      BMM_TLE_EXTS(BMM_TLE_NUM_IO, "Unexpected end");
    case BMM_IO_READ_ERROR:
      return BMM_IO_READ_ERROR;
  }

  switch (num) {
    case BMM_MSG_NUM_PARTS:
      {
        switch (msg_read(&nc->npart, sizeof nc->npart, NULL)) {
          case BMM_IO_READ_EOF:
            BMM_TLE_EXTS(BMM_TLE_NUM_IO, "Unexpected end");
          case BMM_IO_READ_ERROR:
            return BMM_IO_READ_ERROR;
        }

        struct bmm_dem_part parts[BMM_NPART];

        switch (msg_read(parts, sizeof parts, NULL)) {
          case BMM_IO_READ_EOF:
            BMM_TLE_EXTS(BMM_TLE_NUM_IO, "Unexpected end");
          case BMM_IO_READ_ERROR:
            return BMM_IO_READ_ERROR;
        }

        struct bmm_dem_partc partcs[BMM_NPART];

        switch (msg_read(partcs, sizeof partcs, NULL)) {
          case BMM_IO_READ_EOF:
            BMM_TLE_EXTS(BMM_TLE_NUM_IO, "Unexpected end");
          case BMM_IO_READ_ERROR:
            return BMM_IO_READ_ERROR;
        }

        int nerr;

        size_t index[1];

        static float bogus_time = 0.0f;
        bogus_time += 1.0f;

        index[0] = nc->iframe;
        nerr = nc_put_var1_float(nc->ncid, nc->varid_time,
            index, &bogus_time);
        if (nerr != NC_NOERR) {
          BMM_TLE_EXTS(BMM_TLE_NUM_NC, "%s", nc_strerror(nerr));

          return BMM_IO_READ_ERROR;
        }

        float data[BMM_NPART][NDIM];

        // TODO Use `_FillValue`.
        for (size_t ipart = 0; ipart < BMM_NPART; ++ipart)
          for (size_t idim = 0; idim < NDIM; ++idim)
            data[ipart][idim] = ipart >= nc->npart ?
              NAN : (float) parts[ipart].lin.r[idim];

        size_t start[3];
        size_t count[3];

        start[0] = nc->iframe;
        start[1] = 0;
        start[2] = 0;
        count[0] = 1;
        count[1] = BMM_NPART;
        count[2] = NDIM;
        nerr = nc_put_vara_float(nc->ncid, nc->varid_coords,
            start, count, &data[0][0]);
        if (nerr != NC_NOERR) {
          BMM_TLE_EXTS(BMM_TLE_NUM_NC, "%s", nc_strerror(nerr));

          return BMM_IO_READ_ERROR;
        }

        ++nc->iframe;
      }

      break;
    default:
      dynamic_assert(false, "Unsupported message number");
  }

  return BMM_IO_READ_SUCCESS;
}

static bool bmm_nc_run_(struct bmm_nc* const nc) {
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

    switch (bmm_nc_step(nc)) {
      case BMM_IO_READ_ERROR:
        return false;
      case BMM_IO_READ_EOF:
        return true;
    }
  }
}

bool bmm_nc_run(struct bmm_nc* const nc) {
  int const sigs[] = {SIGINT, SIGQUIT, SIGTERM, SIGPIPE};
  if (bmm_sig_register(sigs, nmembof(sigs)) != SIZE_MAX) {
    BMM_TLE_STDS();

    return false;
  }

  if (!bmm_nc_open(nc))
    return false;

  bool const run = bmm_nc_run_(nc);

  if (!bmm_nc_close(nc))
    return false;

  return run;
}

bool bmm_nc_run_with(struct bmm_nc_opts const* const opts) {
  struct bmm_nc nc;
  bmm_nc_def(&nc, opts);

  return bmm_nc_run(&nc);
}
