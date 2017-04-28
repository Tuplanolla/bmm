#ifndef BMM_NC_H
/// Trajectory exporting in NetCDF format.
#define BMM_NC_H

#include <stdbool.h>
#include <stdio.h>

#include "ext.h"

/// This enumeration is for file format conventions.
enum bmm_nc_conv {
  BMM_NC_CONV_AMBER
};

/// This structure contains export options.
struct bmm_nc_opts {
  enum bmm_nc_conv conv;
  char path[BUFSIZ];
  bool i;
  bool r;
  bool q;
  bool v;
  bool f;
};

/// This structure tracks resources.
struct bmm_nc {
  struct bmm_nc_opts opts;
  size_t npart;
  FILE* stream;
  int ncid;
  int id_frame;
  int id_spatial;
  int id_atom;
  int id_cspatial;
  int id_label;
  int varid_spatial;
  int varid_cspatial;
  int varid_time;
  int varid_coords;
  int varid_clens;
  int varid_radii;
  int varid_vels;
  int varid_forces;
};

/// The call `bmm_nc_opts_def(opts)`
/// writes the default export options into `opts`.
/// All messages are stopped by default.
__attribute__ ((__nonnull__))
void bmm_nc_opts_def(struct bmm_nc_opts*);

/// The call `bmm_nc_def(nc)`
/// writes the default export state into `nc`.
__attribute__ ((__nonnull__))
void bmm_nc_def(struct bmm_nc*, struct bmm_nc_opts const*);

/// The call `bmm_nc_step(nc)`
/// processes one incoming message with the export state `nc`.
__attribute__ ((__nonnull__))
enum bmm_io_read bmm_nc_step(struct bmm_nc*);

/// The call `bmm_nc_run(nc)`
/// processes all incoming messages and
/// handles signals with the export state `nc`.
__attribute__ ((__nonnull__))
bool bmm_nc_run(struct bmm_nc*);

/// The call `bmm_nc_run_with(opts)`
/// processes all incoming messages and
/// handles signals with the export options `opts`.
__attribute__ ((__nonnull__))
bool bmm_nc_run_with(struct bmm_nc_opts const*);

#endif
