#include <netcdf.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// clang $(pkg-config --cflags netcdf) -Wall -Wextra -o a.out nc.c $(pkg-config --libs netcdf) && ./a.out && ncdump -h nc.nc

#define URGH(e) ((void) printf("%s\n", nc_strerror(e)))

// Must have `NDIM = 3` for OVITO.
#define NDIM 3
#define NFRAME 12
#define NATOM 64

int main(void) {
  float data[NATOM][NFRAME][NDIM];

  for (int x = 0; x < NATOM; x++)
    for (int y = 0; y < NFRAME; y++)
      for (int z = 0; z < NDIM; z++)
        data[x][y][z] = (float) (rand() % 256) / 256.0f * 10.0f;

  int ncerr;

  int ncid;

  ncerr = nc_create("nc.nc", NC_CLOBBER | NC_64BIT_OFFSET, &ncid);
  if (ncerr != NC_NOERR) {
    URGH(ncerr);

    return 1;
  }

  char const convs[] = "AMBER";
  ncerr = nc_put_att_text(ncid, NC_GLOBAL,
      "Conventions", strlen(convs), convs);
  if (ncerr != NC_NOERR) {
    URGH(ncerr);

    return 1;
  }

  char const convver[] = "1.0";
  ncerr = nc_put_att_text(ncid, NC_GLOBAL,
      "ConventionVersion", strlen(convver), convver);
  if (ncerr != NC_NOERR) {
    URGH(ncerr);

    return 1;
  }

  char const prog[] = "bmm";
  ncerr = nc_put_att_text(ncid, NC_GLOBAL,
      "program", strlen(prog), prog);
  if (ncerr != NC_NOERR) {
    URGH(ncerr);

    return 1;
  }

  char const progver[] = "0.0.0";
  ncerr = nc_put_att_text(ncid, NC_GLOBAL,
      "programVersion", strlen(progver), progver);
  if (ncerr != NC_NOERR) {
    URGH(ncerr);

    return 1;
  }

  int id_frame;
  ncerr = nc_def_dim(ncid, "frame", NC_UNLIMITED, &id_frame);
  if (ncerr != NC_NOERR) {
    URGH(ncerr);

    return 1;
  }

  int id_spatial;
  ncerr = nc_def_dim(ncid, "spatial", NDIM, &id_spatial);
  if (ncerr != NC_NOERR) {
    URGH(ncerr);

    return 1;
  }

  int id_atom;
  ncerr = nc_def_dim(ncid, "atom", NATOM, &id_atom);
  if (ncerr != NC_NOERR) {
    URGH(ncerr);

    return 1;
  }

  int id_cspatial;
  ncerr = nc_def_dim(ncid, "cell_spatial", NDIM, &id_cspatial);
  if (ncerr != NC_NOERR) {
    URGH(ncerr);

    return 1;
  }

  int id_label;
  ncerr = nc_def_dim(ncid, "label", NDIM, &id_label);
  if (ncerr != NC_NOERR) {
    URGH(ncerr);

    return 1;
  }

  int dimids[4];

  int varid_spatial;
  dimids[0] = id_spatial;
  ncerr = nc_def_var(ncid, "spatial", NC_CHAR, 1, dimids, &varid_spatial);
  if (ncerr != NC_NOERR) {
    URGH(ncerr);

    return 1;
  }

  int varid_cspatial;
  dimids[0] = id_cspatial;
  ncerr = nc_def_var(ncid, "cell_spatial", NC_CHAR, 1, dimids, &varid_cspatial);
  if (ncerr != NC_NOERR) {
    URGH(ncerr);

    return 1;
  }

  int varid_time;
  dimids[0] = id_frame;
  ncerr = nc_def_var(ncid, "time", NC_FLOAT, 1, dimids, &varid_time);
  if (ncerr != NC_NOERR) {
    URGH(ncerr);

    return 1;
  }

  char const time[] = "second";
  ncerr = nc_put_att_text(ncid, varid_time, "units", strlen(time), time);
  if (ncerr != NC_NOERR) {
    URGH(ncerr);

    return 1;
  }

  int varid_coords;
  dimids[0] = id_frame;
  dimids[1] = id_atom;
  dimids[2] = id_spatial;
  ncerr = nc_def_var(ncid, "coordinates", NC_FLOAT, 3, dimids, &varid_coords);
  if (ncerr != NC_NOERR) {
    URGH(ncerr);

    return 1;
  }

  char const coords[] = "meter";
  ncerr = nc_put_att_text(ncid, varid_coords, "units", strlen(coords), coords);
  if (ncerr != NC_NOERR) {
    URGH(ncerr);

    return 1;
  }

  int varid_clens;
  dimids[0] = id_frame;
  dimids[1] = id_cspatial;
  ncerr = nc_def_var(ncid, "cell_lengths", NC_FLOAT, 2, dimids, &varid_clens);
  if (ncerr != NC_NOERR) {
    URGH(ncerr);

    return 1;
  }

  char const clens[] = "meter";
  ncerr = nc_put_att_text(ncid, varid_clens, "units", strlen(clens), clens);
  if (ncerr != NC_NOERR) {
    URGH(ncerr);

    return 1;
  }

  int varid_radii;
  dimids[0] = id_frame;
  dimids[1] = id_atom;
  // Must have `name = "radius"` for OVITO.
  ncerr = nc_def_var(ncid, "radius", NC_FLOAT, 2, dimids, &varid_radii);
  if (ncerr != NC_NOERR) {
    URGH(ncerr);

    return 1;
  }

  char const radii[] = "meter";
  ncerr = nc_put_att_text(ncid, varid_radii, "units", strlen(radii), radii);
  if (ncerr != NC_NOERR) {
    URGH(ncerr);

    return 1;
  }

  ncerr = nc_enddef(ncid);
  if (ncerr != NC_NOERR) {
    URGH(ncerr);

    return 1;
  }

  ncerr = nc_put_var_text(ncid, varid_spatial, "xyz");
  if (ncerr != NC_NOERR) {
    URGH(ncerr);

    return 1;
  }

  ncerr = nc_put_var_text(ncid, varid_cspatial, "abc");
  if (ncerr != NC_NOERR) {
    URGH(ncerr);

    return 1;
  }

  size_t start[4];
  size_t count[4];

  start[0] = 0;
  count[0] = NFRAME;
  ncerr = nc_put_vara_float(ncid, varid_time, start, count, &data[0][0][0]);
  if (ncerr != NC_NOERR) {
    URGH(ncerr);

    return 1;
  }

  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  count[0] = NFRAME;
  count[1] = NATOM;
  count[2] = NDIM;
  ncerr = nc_put_vara_float(ncid, varid_coords, start, count, &data[0][0][0]);
  if (ncerr != NC_NOERR) {
    URGH(ncerr);

    return 1;
  }

  float lens[NFRAME][NDIM];
  for (int y = 0; y < NFRAME; y++)
    for (int z = 0; z < NDIM; z++)
      lens[y][z] = 10.0f;

  start[0] = 0;
  start[1] = 0;
  count[0] = NFRAME;
  count[1] = NDIM;
  ncerr = nc_put_vara_float(ncid, varid_clens, start, count, &lens[0][0]);
  if (ncerr != NC_NOERR) {
    URGH(ncerr);

    return 1;
  }

  float rad[NATOM][NFRAME];
  for (int x = 0; x < NATOM; x++)
    for (int y = 0; y < NFRAME; y++)
      rad[x][y] = (float) (rand() % 256) / 256.0f * 0.8f + 0.2f;

  start[0] = 0;
  start[1] = 0;
  count[0] = NFRAME;
  count[1] = NATOM;
  ncerr = nc_put_vara_float(ncid, varid_radii, start, count, &rad[0][0]);
  if (ncerr != NC_NOERR) {
    URGH(ncerr);

    return 1;
  }

  ncerr = nc_close(ncid);
  if (ncerr != NC_NOERR) {
    URGH(ncerr);

    return 1;
  }

  puts("SUCCESS!");

  return 0;
}
