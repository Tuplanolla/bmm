#include <GL/gl.h>
#include <GL/glut.h>
#include <SDL2/SDL.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/time.h>

#include "common.h"
#include "dem.h"
#include "endy.h"
#include "ext.h"
#include "fp.h"
#include "gl.h"
#include "io.h"
#include "msg.h"
#include "sdl.h"
#include "tle.h"

static SDL_Window *window;
static SDL_GLContext glcontext;

extern inline void bmm_sdl_t_to_timeval(struct timeval *, Uint32);

extern inline Uint32 bmm_sdl_t_from_timeval(struct timeval const *);

extern inline Uint32 bmm_sdl_trem(Uint32, Uint32);

void glString(char const *str, int const x, int const y,
    float const *const color, void *font) {
  glColor3fv(color);
  glRasterPos2i(x, y);
  while (*str != '\0')
    glutBitmapCharacter(font, *str++);
}

static enum bmm_io_read msg_read(void *buf, size_t const n,
    __attribute__ ((__unused__)) void *const ptr) {
  return bmm_io_readin(buf, n);
}

enum bmm_io_read bmm_dem_gets_stuff(struct bmm_dem *const dem,
    enum bmm_msg_num const num) {
  switch (num) {
    case BMM_MSG_NUM_ISTEP:
      return msg_read(&dem->time, sizeof dem->time, NULL);
    case BMM_MSG_NUM_OPTS:
      return msg_read(&dem->opts, sizeof dem->opts, NULL);
    case BMM_MSG_NUM_NEIGH:
      switch (msg_read(&dem->cache, sizeof dem->cache, NULL)) {
        case BMM_IO_READ_ERROR:
          return BMM_IO_READ_ERROR;
        case BMM_IO_READ_EOF:
          return BMM_IO_READ_EOF;
      }

      return msg_read(&dem->pair, sizeof dem->pair, NULL);
    case BMM_MSG_NUM_PARTS:
      switch (msg_read(&dem->part.n, sizeof dem->part.n, NULL)) {
        case BMM_IO_READ_ERROR:
          return BMM_IO_READ_ERROR;
        case BMM_IO_READ_EOF:
          return BMM_IO_READ_EOF;
      }

      return msg_read(&dem->part, sizeof dem->part, NULL);
    case BMM_MSG_NUM_EST:
      return msg_read(&dem->est, sizeof dem->est, NULL);
  }

  dynamic_assert(false, "Unsupported message number");
}

enum bmm_io_read bmm_dem_gets(struct bmm_dem *const dem,
    enum bmm_msg_num *const num) {
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

  switch (bmm_msg_num_read(num, msg_read, NULL)) {
    case BMM_IO_READ_ERROR:
      return BMM_IO_READ_ERROR;
    case BMM_IO_READ_EOF:
      return BMM_IO_READ_EOF;
  }

  if (bmm_dem_sniff_size(dem, *num) != spec.msg.size - BMM_MSG_NUMSIZE) {
    BMM_TLE_EXTS(BMM_TLE_NUM_UNKNOWN, "Size mismatch");

    return BMM_IO_READ_ERROR;
  }

  return bmm_dem_gets_stuff(dem, *num);
}

void bmm_sdl_opts_def(struct bmm_sdl_opts *const opts) {
  opts->width = 640;
  opts->height = 480;
  opts->fps = 16;
  opts->ms = 0;
  opts->zoomfac = 1.5;
}

void bmm_sdl_def(struct bmm_sdl *const sdl,
    struct bmm_sdl_opts const *const opts) {
  sdl->opts = *opts;
  sdl->width = (int) opts->width;
  sdl->height = (int) opts->height;
  sdl->qaspect = (double) opts->width / (double) opts->height;
  sdl->qzoom = 1.0;
  sdl->rorigin[0] = 0.0;
  sdl->rorigin[1] = 0.0;
  sdl->fps = opts->fps;
  sdl->itarget = SIZE_MAX;
  sdl->stale = true;
  sdl->active = true;
  sdl->blend = true;
  sdl->diag = true;

  struct bmm_dem_opts defopts;
  bmm_dem_opts_def(&defopts);
  bmm_dem_def(&sdl->dem, &defopts);
}

static Uint32 bmm_sdl_tstep(struct bmm_sdl const *const sdl) {
  return sdl->fps > 1000 ? 1 : (Uint32) (1000 / sdl->fps);
}

// TODO Express this mess in terms of linear algebra.
// The viewport mapping is essentially a homogeneous coordinate transformation.

static void bmm_sdl_proj(struct bmm_sdl const *const sdl,
    double *const xproj, double *const yproj,
    double *const wproj, double *const hproj) {
  double const w = sdl->dem.opts.box.x[0];
  double const h = sdl->dem.opts.box.x[1];
  double const q = sdl->qaspect;
  double const z = sdl->qzoom;
  double const xorigin = sdl->rorigin[0];
  double const yorigin = sdl->rorigin[1];

  bool const pwide = w > h * q;

  double const wzoom = w / z;
  double const hzoom = h / z;
  double const weither = pwide ? wzoom : hzoom * q;
  double const heither = pwide ? wzoom / q : hzoom;

  *wproj = weither;
  *hproj = heither;
  *xproj = xorigin - (weither - wzoom) * 0.5;
  *yproj = yorigin - (heither - hzoom) * 0.5;
}

static void bmm_sdl_zoom(struct bmm_sdl *const sdl,
    double const xscreen, double const yscreen, double const q) {
  double xproj;
  double yproj;
  double wproj;
  double hproj;
  bmm_sdl_proj(sdl, &xproj, &yproj, &wproj, &hproj);

  double const x = bmm_fp_lerp(xscreen,
      0.0, (double) sdl->width, xproj, xproj + wproj);
  double const y = bmm_fp_lerp(yscreen,
      (double) sdl->height, 0.0, yproj, yproj + hproj);

  sdl->rorigin[0] = x;
  sdl->rorigin[1] = y;

  sdl->qzoom *= q;

  bmm_sdl_proj(sdl, &xproj, &yproj, &wproj, &hproj);

  double const x2 = bmm_fp_lerp(xscreen,
      0.0, (double) sdl->width, xproj, xproj + wproj);
  double const y2 = bmm_fp_lerp(yscreen,
      (double) sdl->height, 0.0, yproj, yproj + hproj);

  sdl->rorigin[0] -= x2 - x;
  sdl->rorigin[1] -= y2 - y;
}

static void bmm_sdl_reset(struct bmm_sdl *const sdl) {
  sdl->qzoom = 1.0;
  sdl->rorigin[0] = 0.0;
  sdl->rorigin[1] = 0.0;
}

static void bmm_sdl_move(struct bmm_sdl *const sdl,
    double const xscreen, double const yscreen) {
  double xproj;
  double yproj;
  double wproj;
  double hproj;
  bmm_sdl_proj(sdl, &xproj, &yproj, &wproj, &hproj);

  double const x2 = bmm_fp_lerp(xscreen,
      0.0, (double) sdl->width, xproj, xproj + wproj);
  double const y2 = bmm_fp_lerp(yscreen,
      (double) sdl->height, 0.0, yproj, yproj + hproj);

  sdl->rorigin[0] += x2 - xproj - wproj * 0.5;
  sdl->rorigin[1] += y2 - yproj - hproj * 0.5;
}

void bmm_dem_ijcellx(size_t *const pijcell,
    struct bmm_dem const *const dem, double *y) {
  for (size_t idim = 0; idim < BMM_NDIM; ++idim) {
    size_t const n = dem->opts.cache.ncell[idim];

    size_t const j = n - 1;
    double const r = (double) j;

    double const x = bmm_fp_lerp(y[idim],
        0.0, dem->opts.box.x[idim], 1.0, r);

    if (x < 1.0)
      pijcell[idim] = 0;
    else if (x >= r)
      pijcell[idim] = j;
    else {
      size_t const k = (size_t) x;

      dynamic_assert(k < j, "Invalid truncation");

      pijcell[idim] = k;
    }
  }
}

static void bmm_sdl_draw(struct bmm_sdl const *const sdl) {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

  double xproj;
  double yproj;
  double wproj;
  double hproj;
  bmm_sdl_proj(sdl, &xproj, &yproj, &wproj, &hproj);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(xproj, xproj + wproj, yproj, yproj + hproj, -1.0, 1.0);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  // This would functions as `glClear` with alpha.
  /*
  GLfloat veil[4];
  memcpy(veil, glBlack, sizeof glBlack);
  veil[3] = 0.5f;
  glColor4fv(veil);
  glRectf(xproj, yproj, xproj + wproj, yproj + hproj);
  */

  size_t const ncorner = 8;

  size_t const ncx = sdl->dem.opts.cache.ncell[0] - 2;
  size_t const ncy = sdl->dem.opts.cache.ncell[1] - 2;

  double const w = sdl->dem.opts.box.x[0] / (double) ncx;
  double const h = sdl->dem.opts.box.x[1] / (double) ncy;

  // Bidirectional mappings.
  struct {
    size_t n;
    size_t itgt[BMM_MCONTACT];
  } ind[BMM_MPART];
  for (size_t ipart = 0; ipart < sdl->dem.part.n; ++ipart)
    ind[ipart].n = 0;
  for (size_t ipart = 0; ipart < sdl->dem.part.n; ++ipart)
    for (size_t icont = 0; icont < sdl->dem.pair[BMM_DEM_CT_STRONG].cont.src[ipart].n; ++icont) {
      size_t const jpart = sdl->dem.pair[BMM_DEM_CT_STRONG].cont.src[ipart].itgt[icont];

      size_t const jcont = ind[jpart].n;
      ind[jpart].itgt[jcont] = ipart;
      ++ind[jpart].n;

      // This makes them bidirectional.
      size_t const kcont = ind[ipart].n;
      ind[ipart].itgt[kcont] = jpart;
      ++ind[ipart].n;
    }

  // Particles.
  for (double off = -1; off < 2; ++off) {
    for (size_t ipart = 0; ipart < sdl->dem.part.n; ++ipart) {
      double const x = sdl->dem.part.x[ipart][0];
      double const y = sdl->dem.part.x[ipart][1];
      double const r = sdl->dem.part.r[ipart];
      double const a = sdl->dem.part.phi[ipart];

      double const xoff = x + off * sdl->dem.opts.box.x[0];

      GLfloat blent[4];
      memcpy(blent, sdl->dem.part.role[ipart] == BMM_DEM_ROLE_FREE ?
          glYellow : glWhite, sizeof glBlack);
      blent[3] = 1.0f;
      GLfloat nope[4];
      memcpy(nope, blent, sizeof glBlack);
      nope[3] = 0.0f;
      double t = fabs(xoff / sdl->dem.opts.box.x[0] - 0.5) - 0.5;
      blent[3] = 1.0f - (float) t;
      glColor4fv(blent);

      if (sdl->blend)
        glSkewedAnnulus((float) xoff, (float) y,
            (float) r, (float) (r * 0.25), (float) (r * 0.5),
            (float) a, ncorner);
      else
        glDisk((float) xoff, (float) y, (float) r, ncorner);

      // Cached ghosts.
      if (sdl->blend) {
        double x0[BMM_NDIM];
        double x1[BMM_NDIM];

        (void) memcpy(x0, sdl->dem.part.x[ipart], sizeof x0);
        (void) memcpy(x1, sdl->dem.cache.x[ipart], sizeof x1);

        double dx[BMM_NDIM];
        dx[0] = x0[0] + $(bmm_swrap, double)(x1[0] - x0[0], sdl->dem.opts.box.x[0]);
        dx[1] = x1[1];

        x0[0] += off * sdl->dem.opts.box.x[0];
        dx[0] += off * sdl->dem.opts.box.x[0];

        glBegin(GL_LINES);
        glVertex2dv(x0);
        glVertex2dv(dx);
        glEnd();
      }

      // Contacts.
      if (sdl->blend) {
        memcpy(blent, glRed, sizeof glBlack);
        blent[3] = 1.0f;
        memcpy(nope, blent, sizeof glBlack);
        nope[3] = 0.0f;
        blent[3] = 1.0f - (float) t;
        glColor4fv(blent);

        for (size_t icont = 0; icont < sdl->dem.pair[BMM_DEM_CT_WEAK].cont.src[ipart].n; ++icont) {
          size_t const jpart = sdl->dem.pair[BMM_DEM_CT_WEAK].cont.src[ipart].itgt[icont];

          double x0[BMM_NDIM];
          double x1[BMM_NDIM];

          (void) memcpy(x0, sdl->dem.part.x[ipart], sizeof x0);
          (void) memcpy(x1, sdl->dem.part.x[jpart], sizeof x1);

          double dx[BMM_NDIM];
          dx[0] = x0[0] + $(bmm_swrap, double)(x1[0] - x0[0], sdl->dem.opts.box.x[0]);
          dx[1] = x1[1];

          x0[0] += off * sdl->dem.opts.box.x[0];
          dx[0] += off * sdl->dem.opts.box.x[0];

          glBegin(GL_LINES);
          glVertex2dv(x0);
          glVertex2dv(dx);
          glEnd();
        }

        memcpy(blent, glGreen, sizeof glBlack);
        blent[3] = 1.0f;
        memcpy(nope, blent, sizeof glBlack);
        nope[3] = 0.0f;
        blent[3] = 1.0f - (float) t;
        glColor4fv(blent);

        for (size_t icont = 0; icont < sdl->dem.pair[BMM_DEM_CT_STRONG].cont.src[ipart].n; ++icont) {
          size_t const jpart = sdl->dem.pair[BMM_DEM_CT_STRONG].cont.src[ipart].itgt[icont];

          double x0[BMM_NDIM];
          double x1[BMM_NDIM];

          (void) memcpy(x0, sdl->dem.part.x[ipart], sizeof x0);
          (void) memcpy(x1, sdl->dem.part.x[jpart], sizeof x1);

          double dx[BMM_NDIM];
          dx[0] = x0[0] + $(bmm_swrap, double)(x1[0] - x0[0], sdl->dem.opts.box.x[0]);
          dx[1] = x1[1];

          x0[0] += off * sdl->dem.opts.box.x[0];
          dx[0] += off * sdl->dem.opts.box.x[0];

          glBegin(GL_LINES);
          glVertex2dv(x0);
          glVertex2dv(dx);
          glEnd();
        }
      }

      // Neighbors.
      memcpy(blent, glCyan, sizeof glBlack);
      blent[3] = 1.0f;
      memcpy(nope, blent, sizeof glBlack);
      nope[3] = 0.0f;
      blent[3] = 1.0f - (float) t;
      glColor4fv(blent);

      // Focus.
      if (ipart != sdl->itarget)
        continue;

      glBegin(GL_LINES);
      for (size_t ineigh = 0; ineigh < sdl->dem.cache.neigh[ipart].n; ++ineigh) {
        size_t const jpart = sdl->dem.cache.neigh[ipart].i[ineigh];

        double x0[BMM_NDIM];
        double x1[BMM_NDIM];

        (void) memcpy(x0, sdl->dem.part.x[ipart], sizeof x0);
        (void) memcpy(x1, sdl->dem.part.x[jpart], sizeof x1);

        double dx[BMM_NDIM];
        dx[0] = x0[0] + $(bmm_swrap, double)(x1[0] - x0[0], sdl->dem.opts.box.x[0]);
        dx[1] = x1[1];

        x0[0] += off * sdl->dem.opts.box.x[0];
        dx[0] += off * sdl->dem.opts.box.x[0];

        glVertex2dv(x0);
        glVertex2dv(dx);
      }
      glEnd();
    }

    // Beams.
    if (!sdl->blend) {
      glBegin(GL_QUADS);
      for (size_t ipart = 0; ipart < sdl->dem.part.n; ++ipart)
        for (size_t icont = 0; icont < sdl->dem.pair[BMM_DEM_CT_STRONG].cont.src[ipart].n; ++icont) {
          size_t const jpart = sdl->dem.pair[BMM_DEM_CT_STRONG].cont.src[ipart].itgt[icont];

          double xdiffji[BMM_NDIM];
          bmm_geom2d_cpdiff(xdiffji, sdl->dem.part.x[ipart], sdl->dem.part.x[jpart],
              sdl->dem.opts.box.x, sdl->dem.opts.box.per);

          double const d2 = bmm_geom2d_norm2(xdiffji);
          double const ri = sdl->dem.part.r[ipart];
          double const rj = sdl->dem.part.r[jpart];
          double const d = sqrt(d2);

          double xnormji[BMM_NDIM];
          bmm_geom2d_scale(xnormji, xdiffji, 1.0 / d);

          double xtangji[BMM_NDIM];
          bmm_geom2d_rperp(xtangji, xnormji);

          double x0[BMM_NDIM];
          double x1[BMM_NDIM];
          double x2[BMM_NDIM];
          double x3[BMM_NDIM];
          bmm_geom2d_scale(x0, xtangji, ri);
          bmm_geom2d_scale(x1, x0, -1.0);
          bmm_geom2d_scale(x2, xtangji, rj);
          bmm_geom2d_scale(x3, x2, -1.0);

          bmm_geom2d_addto(x0, sdl->dem.part.x[ipart]);
          bmm_geom2d_addto(x1, sdl->dem.part.x[ipart]);
          bmm_geom2d_addto(x2, sdl->dem.part.x[jpart]);
          bmm_geom2d_addto(x3, sdl->dem.part.x[jpart]);

          double dx[BMM_NDIM];
          dx[0] = x0[0] + $(bmm_swrap, double)(x1[0] - x0[0], sdl->dem.opts.box.x[0]);
          dx[1] = x1[1];

          double dy[BMM_NDIM];
          dy[0] = x0[0] + $(bmm_swrap, double)(x2[0] - x0[0], sdl->dem.opts.box.x[0]);
          dy[1] = x2[1];

          double dz[BMM_NDIM];
          dz[0] = x0[0] + $(bmm_swrap, double)(x3[0] - x0[0], sdl->dem.opts.box.x[0]);
          dz[1] = x3[1];

          x0[0] += off * sdl->dem.opts.box.x[0];
          dx[0] += off * sdl->dem.opts.box.x[0];
          dy[0] += off * sdl->dem.opts.box.x[0];
          dz[0] += off * sdl->dem.opts.box.x[0];

          GLfloat blent[4];
          memcpy(blent, sdl->dem.part.role[ipart] == BMM_DEM_ROLE_FREE ?
              glYellow : glWhite, sizeof glBlack);
          blent[3] = 0.5f;

          GLfloat blunt[4];
          memcpy(blunt, sdl->dem.part.role[jpart] == BMM_DEM_ROLE_FREE ?
              glYellow : glWhite, sizeof glBlack);
          blunt[3] = 0.5f;

          glColor4fv(blent);
          glVertex2dv(x0);
          glColor4fv(blunt);
          glVertex2dv(dy);
          glVertex2dv(dz);
          glColor4fv(blent);
          glVertex2dv(dx);
        }
      glEnd();
    }

    // Connected components (three deep).
    if (false && sdl->blend) {
      glBegin(GL_TRIANGLES);
      for (size_t ipart = 0; ipart < sdl->dem.part.n; ++ipart)
        for (size_t icont = 0; icont < ind[ipart].n; ++icont) {
          size_t const jpart = ind[ipart].itgt[icont];

          for (size_t jcont = 0; jcont < ind[jpart].n; ++jcont) {
            size_t const kpart = ind[jpart].itgt[jcont];

            for (size_t kcont = 0; kcont < ind[kpart].n; ++kcont) {
              size_t const lpart = ind[kpart].itgt[kcont];

              if (lpart == ipart &&
                  sdl->dem.part.x[ipart][0] < sdl->dem.part.x[jpart][0] &&
                  sdl->dem.part.x[ipart][0] < sdl->dem.part.x[kpart][0]) {

                double x0[BMM_NDIM];
                double x1[BMM_NDIM];
                double x2[BMM_NDIM];

                (void) memcpy(x0, sdl->dem.part.x[ipart], sizeof x0);
                (void) memcpy(x1, sdl->dem.part.x[jpart], sizeof x1);
                (void) memcpy(x2, sdl->dem.part.x[kpart], sizeof x2);

                double dx[BMM_NDIM];
                dx[0] = x0[0] + $(bmm_swrap, double)(x1[0] - x0[0], sdl->dem.opts.box.x[0]);
                dx[1] = x1[1];

                double dy[BMM_NDIM];
                dy[0] = x0[0] + $(bmm_swrap, double)(x2[0] - x0[0], sdl->dem.opts.box.x[0]);
                dy[1] = x2[1];

                x0[0] += off * sdl->dem.opts.box.x[0];
                dx[0] += off * sdl->dem.opts.box.x[0];
                dy[0] += off * sdl->dem.opts.box.x[0];

                {
                  double const x = sdl->dem.part.x[ipart][0];
                  double const y = sdl->dem.part.x[ipart][1];

                  double const xoff = x + off * sdl->dem.opts.box.x[0];

                  GLfloat blent[4];
                  memcpy(blent, sdl->dem.part.role[ipart] == BMM_DEM_ROLE_FREE ?
                      glYellow : glWhite, sizeof glBlack);
                  blent[3] = 1.0f;
                  GLfloat nope[4];
                  memcpy(nope, blent, sizeof glBlack);
                  nope[3] = 0.0f;
                  double t = fabs(xoff / sdl->dem.opts.box.x[0] - 0.5) - 0.5;
                  blent[3] = 1.0f - (float) t;
                  glColor4fv(blent);

                  glVertex2dv(x0);
                  glVertex2dv(dx);
                  glVertex2dv(dy);
                }
              }
            }
          }
        }
      glEnd();
    }
  }

  // Staleness indicator.
  if (sdl->diag) {
    glColor3fv(sdl->stale ? glRed : glGreen);
    glDisk((float) (0.05 * sdl->dem.opts.box.x[0]),
        (float) (0.05 * sdl->dem.opts.box.x[1]),
        (float) (0.025 * sdl->dem.opts.box.x[0]), ncorner);
  }

  // Cell boxes.
  if (sdl->blend) {
    glColor3fv(glCyan);
    for (size_t icellx = 0; icellx < ncx; ++icellx)
      for (size_t icelly = 0; icelly < ncy; ++icelly) {
        double const x = (double) icellx * w;
        double const y = (double) icelly * h;

        glRectWire((float) x, (float) y, (float) w, (float) h);
      }
  }

  // Bounding box.
  glColor3fv(glWhite);
  glRectWire(0.0f, 0.0f, (float) sdl->dem.opts.box.x[0], (float) sdl->dem.opts.box.x[1]);

  // Diagnostic text.
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0.0, sdl->width, sdl->height, 0.0, -1.0, 1.0);

  struct bmm_dem const *const dem = &sdl->dem;

  double const eambdis = dem->est.eambdis;
  double const epotext = dem->est.epotext;
  double const eklin = dem->est.eklin;
  double const ekrot = dem->est.ekrot;
  double const ewcont = dem->est.ewcont;
  double const escont = dem->est.escont;
  double const edrivnorm = dem->est.edrivnorm;
  double const edrivtang = dem->est.edrivtang;
  double const eyieldis = dem->est.eyieldis;
  double const ewcontdis = dem->est.ewcontdis;
  double const escontdis = dem->est.escontdis;
  double const pos = eambdis + epotext + eklin + ekrot + ewcont + escont
    + eyieldis + ewcontdis + escontdis;
  double const neg = edrivnorm + edrivtang;
  double const eee = pos - neg;

  // TODO These should come via messages.
  if (sdl->diag) {
    char strbuf[BUFSIZ];
    int ioff = 1;
    (void) snprintf(strbuf, sizeof strbuf, "f (target vfps) = %u (%u)",
        sdl->fps, bmm_sdl_tstep(sdl));
    glString(strbuf, 8, 8 + 15 * ioff++, glWhite, GLUT_BITMAP_9_BY_15);
    (void) snprintf(strbuf, sizeof strbuf, "t (now) = %g (%zu)",
        sdl->dem.time.t, sdl->dem.time.istep);
    glString(strbuf, 8, 8 + 15 * ioff++, glWhite, GLUT_BITMAP_9_BY_15);
    (void) snprintf(strbuf, sizeof strbuf, "t (prev. cache refresh) = %g",
        sdl->dem.cache.tprev);
    glString(strbuf, 8, 8 + 15 * ioff++, glWhite, GLUT_BITMAP_9_BY_15);
    (void) snprintf(strbuf, sizeof strbuf, "n (number of particles) = %zu",
        sdl->dem.part.n);
    glString(strbuf, 8, 8 + 15 * ioff++, glWhite, GLUT_BITMAP_9_BY_15);
    (void) snprintf(strbuf, sizeof strbuf, "W (work in) = %g",
        neg);
    glString(strbuf, 8, 8 + 15 * ioff++, glWhite, GLUT_BITMAP_9_BY_15);
    (void) snprintf(strbuf, sizeof strbuf, "E (energy balance) = %g",
        eee);
    glString(strbuf, 8, 8 + 15 * ioff++, glWhite, GLUT_BITMAP_9_BY_15);
    (void) snprintf(strbuf, sizeof strbuf, "mu (eff. friction factor) = %g",
        sdl->dem.est.mueff);
    glString(strbuf, 8, 8 + 15 * ioff++, glWhite, GLUT_BITMAP_9_BY_15);
    (void) snprintf(strbuf, sizeof strbuf, "F (force feedback) = (%g, %g)",
        sdl->dem.est.fback[0], sdl->dem.est.fback[1]);
    glString(strbuf, 8, 8 + 15 * ioff++, glWhite, GLUT_BITMAP_9_BY_15);
    (void) snprintf(strbuf, sizeof strbuf, "v (driving velocity) = (%g, %g)",
        sdl->dem.est.vdriv[0], sdl->dem.est.vdriv[1]);
    glString(strbuf, 8, 8 + 15 * ioff++, glWhite, GLUT_BITMAP_9_BY_15);
  }

  SDL_GL_SwapWindow(window);
}

static bool bmm_sdl_video(struct bmm_sdl *const sdl,
    int const width, int const height) {
  if (window == NULL) {
    Uint32 const flags = SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE;
    window = SDL_CreateWindow("BMM",
        SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
        width, height, flags);
    if (window == NULL) {
      BMM_TLE_EXTS(BMM_TLE_NUM_SDL, "SDL error: %s", SDL_GetError());

      return false;
    }
  }

  glcontext = SDL_GL_CreateContext(window);

  if (sdl->blend)
    glEnable(GL_BLEND);
  // glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE);
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);

  glClearColor3fv(glBlack);

  glViewport(0, 0, width, height);

  sdl->width = width;
  sdl->height = height;
  sdl->qaspect = (double) width / (double) height;

  return true;
}

// TODO Move heresy.

static bool heresy(struct bmm_sdl const *const sdl) {
  FILE *const stream = fopen("heresy.data", "w");
  if (stream == NULL) {
    BMM_TLE_STDS();

    return false;
  }

  for (size_t ipart = 0; ipart < sdl->dem.part.n; ++ipart)
    if (fprintf(stream, "%zu %g %g %g\n",
          ipart,
          sdl->dem.part.x[ipart][0],
          sdl->dem.part.x[ipart][1],
          sdl->dem.part.r[ipart]) < 0) {
      BMM_TLE_STDS();

      break;
    }

  if (fclose(stream) != 0) {
    BMM_TLE_STDS();

    return false;
  }

  return true;
}

static bool more_heresy(struct bmm_sdl const *const sdl) {
  FILE *const stream = fopen("more-heresy.data", "w");
  if (stream == NULL) {
    BMM_TLE_STDS();

    return false;
  }

  for (size_t ipart = 0; ipart < sdl->dem.part.n; ++ipart) {
    for (size_t icont = 0; icont < sdl->dem.pair[BMM_DEM_CT_STRONG].cont.src[ipart].n; ++icont) {
      size_t const jpart = sdl->dem.pair[BMM_DEM_CT_STRONG].cont.src[ipart].itgt[icont];

      if (fprintf(stream, "%zu %g %g\n",
            (size_t) 0,
            sdl->dem.part.x[ipart][0],
            sdl->dem.part.x[ipart][1]) < 0) {
        BMM_TLE_STDS();

        break; // out
      }

      if (fprintf(stream, "%zu %g %g\n",
            (size_t) 1,
            sdl->dem.part.x[jpart][0],
            sdl->dem.part.x[jpart][1]) < 0) {
        BMM_TLE_STDS();

        break; // out
      }

      if (fprintf(stream, "\n") < 0) {
        BMM_TLE_STDS();

        break; // out
      }
    }
  }

  if (fclose(stream) != 0) {
    BMM_TLE_STDS();

    return false;
  }

  return true;
}

// TODO This is a shit file format.

static bool serious_heresy(struct bmm_sdl const *const sdl) {
  FILE *const stream = fopen("heresy.xyz", "w");
  if (stream == NULL) {
    BMM_TLE_STDS();

    return false;
  }

  if (fprintf(stream, "%zu\n.\n", sdl->dem.part.n) < 0) {
    BMM_TLE_STDS();

    return false;
  }

  for (size_t ipart = 0; ipart < sdl->dem.part.n; ++ipart)
    if (fprintf(stream, "S %g %g 0.0 %g\n",
          sdl->dem.part.x[ipart][0],
          sdl->dem.part.x[ipart][1],
          sdl->dem.part.r[ipart]) < 0) {
      BMM_TLE_STDS();

      break;
    }

  if (fclose(stream) != 0) {
    BMM_TLE_STDS();

    return false;
  }

  return true;
}

static bool bmm_sdl_work(struct bmm_sdl *const sdl) {
  if (!bmm_sdl_video(sdl, sdl->width, sdl->height))
    return false;

  Uint32 tnow = SDL_GetTicks();
  Uint32 tnext = tnow + bmm_sdl_tstep(sdl);
  Uint32 trem = bmm_sdl_trem(tnow, tnext);

  for ever {
    // Handle events just before drawing for maximal responsiveness.
    SDL_Event event;
    while (SDL_PollEvent(&event))
      switch (event.type) {
        case SDL_QUIT:
          return true;
        case SDL_WINDOWEVENT:
          switch (event.window.event) {
            case SDL_WINDOWEVENT_RESIZED:
              if (!bmm_sdl_video(sdl, event.window.data1, event.window.data2))
                return false;
              break;
          }
        case SDL_KEYDOWN:
          switch (event.key.keysym.sym) {
            case SDLK_ESCAPE:
            case SDLK_q:
              return heresy(sdl) && more_heresy(sdl) && serious_heresy(sdl);
              // return true;
            case SDLK_SPACE:
              sdl->active = !sdl->active;
              break;
            case SDLK_PAGEDOWN:
              if (sdl->fps % 2 == 0)
                sdl->fps /= 2;
              break;
            case SDLK_PAGEUP:
              sdl->fps *= 2;
              break;
            case SDLK_MINUS:
            case SDLK_KP_MINUS:
              bmm_sdl_zoom(sdl,
                  (double) sdl->width * 0.5, (double) sdl->height * 0.5,
                  1.0 / sdl->opts.zoomfac);
              break;
            case SDLK_PLUS:
            case SDLK_KP_PLUS:
              bmm_sdl_zoom(sdl,
                  (double) sdl->width * 0.5, (double) sdl->height * 0.5,
                  sdl->opts.zoomfac);
              break;
            case SDLK_LEFT:
              bmm_sdl_move(sdl,
                  (double) sdl->width * 0.25, (double) sdl->height * 0.5);
              break;
            case SDLK_RIGHT:
              bmm_sdl_move(sdl,
                  (double) sdl->width * 0.75, (double) sdl->height * 0.5);
              break;
            case SDLK_DOWN:
              bmm_sdl_move(sdl,
                  (double) sdl->width * 0.5, (double) sdl->height * 0.75);
              break;
            case SDLK_UP:
              bmm_sdl_move(sdl,
                  (double) sdl->width * 0.5, (double) sdl->height * 0.25);
              break;
            case SDLK_0:
            case SDLK_KP_0:
              bmm_sdl_reset(sdl);
              break;
            case SDLK_1:
              sdl->itarget = SIZE_MAX;
              break;
            case SDLK_2:
              sdl->itarget = 0;
              break;
            case SDLK_3:
              sdl->itarget = (size_t) rand() % sdl->dem.part.n;
              break;
            case SDLK_b:
              sdl->blend = !sdl->blend;
              if (sdl->blend)
                glEnable(GL_BLEND);
              else
                glDisable(GL_BLEND);
              break;
            case SDLK_d:
              sdl->diag = !sdl->diag;
              break;
          }
          break;
        case SDL_MOUSEBUTTONDOWN:
          switch (event.button.button) {
            case SDL_BUTTON_LEFT:
              bmm_sdl_move(sdl,
                  (double) event.button.x, (double) event.button.y);
              break;
            case SDL_BUTTON_MIDDLE:
              bmm_sdl_reset(sdl);
              break;
          }
          break;
        case SDL_MOUSEWHEEL:
          {
            int x;
            int y;
            SDL_GetMouseState(&x, &y);

            if (event.wheel.y < 0) {
              bmm_sdl_zoom(sdl, (double) x, (double) y,
                  1.0 / sdl->opts.zoomfac);
            } else if (event.wheel.y > 0) {
              bmm_sdl_zoom(sdl, (double) x, (double) y,
                  sdl->opts.zoomfac);
            }
          }
      }

    // Draw before state transformations to account for the initial state.
    bmm_sdl_draw(sdl);

    // Recompute the remaining time after event handling and drawing
    // in case they take a while to complete.
    tnow = SDL_GetTicks();
    trem = bmm_sdl_trem(tnow, tnext);

    // TODO Think about a better way to express this thing.
    // TODO Currently this does staleness analysis wrong.

    if (sdl->active) {
      // Use the remaining time to wait for input.
      // If there is no time left,
      // try to read just one message to prevent congestion.
      struct timeval timeout;
again:
      bmm_sdl_t_to_timeval(&timeout, trem);
      enum bmm_msg_num num;
      switch (bmm_io_waitin(&timeout)) {
        case BMM_IO_WAIT_ERROR:
          return false;
        case BMM_IO_WAIT_READY:
          (void) bmm_dem_gets(&sdl->dem, &num);

          if (num == BMM_MSG_NUM_PARTS)
            sdl->stale = false;
          else {
            trem = bmm_sdl_t_from_timeval(&timeout);
            if (trem > 0)
              goto again;
          }

          break;
        case BMM_IO_WAIT_TIMEOUT:
          sdl->stale = true;
      }
    }

    // Recompute the remaining time again after waiting for input
    // in case it arrives early or takes a while to process.
    tnow = SDL_GetTicks();
    trem = bmm_sdl_trem(tnow, tnext);

    // Sleep the remaining time and allow the tick counter to wrap.
    SDL_Delay(trem);
    tnext += bmm_sdl_tstep(sdl);

    // TODO Not here.
    {
      static int skip = 1;
      if (skip == 0) {
        SDL_Event e;
        e.type = SDL_KEYDOWN;
        e.key.keysym.sym = SDLK_MINUS;
        SDL_PushEvent(&e);
      }
      if (skip >= 0)
        --skip;
    }
  }
}

bool bmm_sdl_run(struct bmm_sdl *const sdl) {
  if (SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 8) == -1 ||
      SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE, 8) == -1 ||
      SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE, 8) == -1 ||
      SDL_GL_SetAttribute(SDL_GL_ALPHA_SIZE, 8) == -1 ||
      SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 8) == -1 ||
      SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1) == -1) {
    BMM_TLE_EXTS(BMM_TLE_NUM_SDL, "SDL error: %s", SDL_GetError());

    return false;
  }

  if (sdl->opts.ms > 0)
    if (SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS, 1) == -1 ||
        SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES,
          (int) sdl->opts.ms) == -1) {
      BMM_TLE_EXTS(BMM_TLE_NUM_SDL, "SDL error: %s", SDL_GetError());

      return false;
    }

  // TODO This is a bad idea and should not even work.
  int argc = 1;
  char arg[] = "";
  char *argv[] = {arg};
  glutInit(&argc, argv);

  if (!bmm_sdl_work(sdl))
    return false;

  return true;
}

static bool bmm_sdl_run_sdl(struct bmm_sdl_opts const *const opts) {
  struct bmm_sdl *const sdl = malloc(sizeof *sdl);
  if (sdl == NULL) {
    BMM_TLE_STDS();

    return false;
  }

  bmm_sdl_def(sdl, opts);
  bool const result = bmm_sdl_run(sdl);

  free(sdl);

  return result;
}

bool bmm_sdl_run_with(struct bmm_sdl_opts const *const opts) {
  if (SDL_Init(SDL_INIT_VIDEO) == -1) {
    BMM_TLE_EXTS(BMM_TLE_NUM_SDL, "SDL error: %s", SDL_GetError());

    return false;
  }

  bool const result = bmm_sdl_run_sdl(opts);

  SDL_Quit();

  return result;
}
