#include <GL/gl.h>
#include <GL/glut.h>
#include <SDL2/SDL.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/time.h>

#include "dem.h"
#include "endy.h"
#include "ext.h"
#include "fp.h"
#include "gl.h"
#include "io.h"
#include "msg.h"
#include "sdl.h"
#include "size.h"
#include "tle.h"

static SDL_Window* window;
static SDL_GLContext glcontext;

extern inline void bmm_sdl_t_to_timeval(struct timeval*, Uint32);

extern inline Uint32 bmm_sdl_t_from_timeval(struct timeval const*);

extern inline Uint32 bmm_sdl_trem(Uint32, Uint32);

void glString(char const* str, int const x, int const y,
    float const* const color, void* font) {
  glColor4fv(color);
  glRasterPos2i(x, y);
  while (*str != '\0')
    glutBitmapCharacter(font, *str++);
}

static enum bmm_io_read msg_read(void* buf, size_t const n,
    __attribute__ ((__unused__)) void* const ptr) {
  return bmm_io_readin(buf, n);
}

enum bmm_io_read bmm_dem_gets_stuff(struct bmm_dem* const dem,
    enum bmm_msg_num const num) {
  switch (num) {
    case BMM_MSG_NUM_NPART:
      return msg_read(&dem->buf.npart, sizeof dem->buf.npart, NULL);
    case BMM_MSG_NUM_EKINE:
      switch (msg_read(&dem->istep, sizeof dem->istep, NULL)) {
        case BMM_IO_READ_ERROR:
          return BMM_IO_READ_ERROR;
        case BMM_IO_READ_EOF:
          return BMM_IO_READ_EOF;
      }

      return msg_read(&dem->est, sizeof dem->est, NULL);
    case BMM_MSG_NUM_NEIGH:
      switch (msg_read(&dem->buf.neigh, sizeof dem->buf.neigh, NULL)) {
        case BMM_IO_READ_ERROR:
          return BMM_IO_READ_ERROR;
        case BMM_IO_READ_EOF:
          return BMM_IO_READ_EOF;
      }

      return msg_read(&dem->buf.links, sizeof dem->buf.links, NULL);
    case BMM_MSG_NUM_PARTS:
      switch (msg_read(&dem->istep, sizeof dem->istep, NULL)) {
        case BMM_IO_READ_ERROR:
          return BMM_IO_READ_ERROR;
        case BMM_IO_READ_EOF:
          return BMM_IO_READ_EOF;
      }

      switch (msg_read(&dem->buf.parts, sizeof dem->buf.parts, NULL)) {
        case BMM_IO_READ_ERROR:
          return BMM_IO_READ_ERROR;
        case BMM_IO_READ_EOF:
          return BMM_IO_READ_EOF;
      }

      return msg_read(&dem->buf.partcs, sizeof dem->buf.partcs, NULL);
  }

  dynamic_assert(false, "Unsupported message num");
}

enum bmm_io_read bmm_dem_gets(struct bmm_dem* const dem,
    enum bmm_msg_num* const num) {
  struct bmm_msg_spec spec;
  switch (bmm_msg_spec_read(&spec, msg_read, NULL)) {
    case BMM_IO_READ_ERROR:
      return BMM_IO_READ_ERROR;
    case BMM_IO_READ_EOF:
      return BMM_IO_READ_EOF;
  }

  if (spec.endy != bmm_endy_get()) {
    BMM_TLE_EXTS(BMM_TLE_UNIMPL, "Unsupported endianness");

    return BMM_IO_READ_ERROR;
  }

  if (spec.tag != BMM_MSG_TAG_SP) {
    BMM_TLE_EXTS(BMM_TLE_UNIMPL, "Unsupported tag");

    return BMM_IO_READ_ERROR;
  }

  switch (bmm_msg_num_read(num, msg_read, NULL)) {
    case BMM_IO_READ_ERROR:
      return BMM_IO_READ_ERROR;
    case BMM_IO_READ_EOF:
      return BMM_IO_READ_EOF;
  }

  if (bmm_dem_sniff_size(dem, *num) != spec.msg.size - BMM_MSG_NUMSIZE) {
    BMM_TLE_EXTS(BMM_TLE_UNKNOWN, "Size mismatch");

    return BMM_IO_READ_ERROR;
  }

  return bmm_dem_gets_stuff(dem, *num);
}

void bmm_sdl_opts_def(struct bmm_sdl_opts* const opts) {
  opts->width = 640;
  opts->height = 480;
  opts->fps = 16;
  opts->ms = 0;
  opts->zoomfac = 1.5;
}

void bmm_sdl_def(struct bmm_sdl* const sdl,
    struct bmm_sdl_opts const* const opts) {
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

  struct bmm_dem_opts defopts;
  bmm_dem_opts_def(&defopts);
  bmm_dem_def(&sdl->dem, &defopts);
}

static Uint32 bmm_sdl_tstep(struct bmm_sdl const* const sdl) {
  return sdl->fps > 1000 ? 1 : (Uint32) (1000 / sdl->fps);
}

// TODO Express this mess in terms of linear algebra.
// The viewport mapping is essentially a homogeneous coordinate transformation.

static void bmm_sdl_proj(struct bmm_sdl const* const sdl,
    double* const xproj, double* const yproj,
    double* const wproj, double* const hproj) {
  double const w = sdl->dem.rext[0];
  double const h = sdl->dem.rext[1];
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

static void bmm_sdl_zoom(struct bmm_sdl* const sdl,
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

static void bmm_sdl_reset(struct bmm_sdl* const sdl) {
  sdl->qzoom = 1.0;
  sdl->rorigin[0] = 0.0;
  sdl->rorigin[1] = 0.0;
}

static void bmm_sdl_move(struct bmm_sdl* const sdl,
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

static void bmm_sdl_draw(struct bmm_sdl const* const sdl) {
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

  size_t const ncorner = 8;

  double const w = sdl->dem.rext[0] / (double) sdl->dem.opts.ncell[0];
  double const h = sdl->dem.rext[1] / (double) sdl->dem.opts.ncell[1];

  // Particles.
  for (size_t ipart = 0; ipart < sdl->dem.buf.npart; ++ipart) {
    float const x = (float) sdl->dem.buf.parts[ipart].lin.r[0];
    float const y = (float) sdl->dem.buf.parts[ipart].lin.r[1];
    float const r = (float) sdl->dem.buf.partcs[ipart].rrad;
    float const a = (float) sdl->dem.buf.parts[ipart].ang.alpha;

    glColor4fv(sdl->dem.buf.partcs[ipart].free ? glYellow : glWhite);

    glSkewedAnnulus(x, y, r, r * 0.25f, r * 0.5f, a, ncorner);
  }

  // Focus.
  if (sdl->itarget != SIZE_MAX) {
    glColor4fv(glBlue);
    size_t const ipart = sdl->itarget;

    // Neighbor markers.
    glBegin(GL_LINES);
    for (size_t ineigh = 0; ineigh < bmm_dem_sizey(&sdl->dem.buf.neigh.neighs[ipart]); ++ineigh) {
      size_t const jpart = bmm_dem_gety(&sdl->dem.buf.neigh.neighs[ipart], ineigh);

      bool p = true;
      for (size_t idim = 0; idim < 2; ++idim)
        p = p &&
          fabs(sdl->dem.buf.parts[ipart].lin.r[idim] - sdl->dem.buf.parts[jpart].lin.r[idim]) <
          sdl->dem.rext[idim] / 2.0;

      if (p) {
        glVertex2dv(sdl->dem.buf.parts[ipart].lin.r);
        glVertex2dv(sdl->dem.buf.parts[jpart].lin.r);
      }
    }
    glEnd();
  }

  // Structural links.
  glColor4fv(glMagenta);
  for (size_t ipart = 0; ipart < sdl->dem.buf.npart; ++ipart) {
    glBegin(GL_LINES);
    for (size_t ineigh = 0; ineigh < bmm_dem_sizel(&sdl->dem.buf.links[ipart]); ++ineigh) {
      size_t const jpart = bmm_dem_getl(&sdl->dem.buf.links[ipart], ineigh);

      bool p = true;
      for (size_t idim = 0; idim < 2; ++idim)
        p = p &&
          fabs(sdl->dem.buf.parts[ipart].lin.r[idim] - sdl->dem.buf.parts[jpart].lin.r[idim]) <
          sdl->dem.rext[idim] / 2.0;

      if (p) {
        glVertex2dv(sdl->dem.buf.parts[ipart].lin.r);
        glVertex2dv(sdl->dem.buf.parts[jpart].lin.r);
      }
    }
    glEnd();
  }

  // Staleness indicator.
  glColor4fv(sdl->stale ? glRed : glGreen);
  glDisk(0.05f, 0.05f, 0.025f, ncorner);

  // Cell boxes.
  glColor4fv(glCyan);
  for (size_t icellx = 0; icellx < sdl->dem.opts.ncell[0]; ++icellx)
    for (size_t icelly = 0; icelly < sdl->dem.opts.ncell[1]; ++icelly) {
      double const x = (double) icellx * w;
      double const y = (double) icelly * h;

      glRectWire((float) x, (float) y, (float) w, (float) h);
    }

  // Bounding box.
  glColor4fv(glWhite);
  glRectWire(0.0f, 0.0f, (float) sdl->dem.rext[0], (float) sdl->dem.rext[1]);

  // Diagnostic text.
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0.0, sdl->width, sdl->height, 0.0, -1.0, 1.0);

  // TODO These should come via messages.
  char strbuf[BUFSIZ];
  int ioff = 1;
  (void) snprintf(strbuf, sizeof strbuf, "f (target vfps) = %u (%u)",
      sdl->fps, bmm_sdl_tstep(sdl));
  glString(strbuf, 8, 8 + 15 * ioff++, glWhite, GLUT_BITMAP_9_BY_15);
  (void) snprintf(strbuf, sizeof strbuf, "K (kinetic energy) = %g",
      bmm_dem_ekinetic(&sdl->dem));
  glString(strbuf, 8, 8 + 15 * ioff++, glWhite, GLUT_BITMAP_9_BY_15);
  (void) snprintf(strbuf, sizeof strbuf, "p (total vector momentum) = %g",
      bmm_dem_pvector(&sdl->dem));
  glString(strbuf, 8, 8 + 15 * ioff++, glWhite, GLUT_BITMAP_9_BY_15);
  (void) snprintf(strbuf, sizeof strbuf, "p (total scalar momentum) = %g",
      bmm_dem_pscalar(&sdl->dem));
  glString(strbuf, 8, 8 + 15 * ioff++, glWhite, GLUT_BITMAP_9_BY_15);
  (void) snprintf(strbuf, sizeof strbuf, "t (now) = %g",
      sdl->dem.istep * sdl->dem.opts.tstep);
  glString(strbuf, 8, 8 + 15 * ioff++, glWhite, GLUT_BITMAP_9_BY_15);
  (void) snprintf(strbuf, sizeof strbuf, "t (next sched. update) = %g",
      sdl->dem.buf.neigh.tnext);
  glString(strbuf, 8, 8 + 15 * ioff++, glWhite, GLUT_BITMAP_9_BY_15);
  (void) snprintf(strbuf, sizeof strbuf, "e (coeff. of restit.) = %g",
      bmm_dem_cor(&sdl->dem));
  glString(strbuf, 8, 8 + 15 * ioff++, glWhite, GLUT_BITMAP_9_BY_15);
  (void) snprintf(strbuf, sizeof strbuf, "n (number of particles) = %zu",
      sdl->dem.buf.npart);
  glString(strbuf, 8, 8 + 15 * ioff++, glWhite, GLUT_BITMAP_9_BY_15);

  SDL_GL_SwapWindow(window);
}

static bool bmm_sdl_video(struct bmm_sdl* const sdl,
    int const width, int const height) {
  if (window == NULL) {
    Uint32 const flags = SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE;
    window = SDL_CreateWindow("BMM",
        SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
        width, height, flags);
    if (window == NULL) {
      BMM_TLE_EXTS(BMM_TLE_SDL, "SDL error: %s", SDL_GetError());

      return false;
    }
  }

  glcontext = SDL_GL_CreateContext(window);

  glViewport(0, 0, width, height);

  sdl->width = width;
  sdl->height = height;
  sdl->qaspect = (double) width / (double) height;

  return true;
}

// TODO Move heresy.

static bool heresy(struct bmm_sdl const* const sdl) {
  FILE* const stream = fopen("heresy.data", "w");
  if (stream == NULL) {
    BMM_TLE_STDS();

    return false;
  }

  for (size_t ipart = 0; ipart < sdl->dem.buf.npart; ++ipart)
    if (fprintf(stream, "%zu %g %g %g\n",
          ipart,
          sdl->dem.buf.parts[ipart].lin.r[0],
          sdl->dem.buf.parts[ipart].lin.r[1],
          sdl->dem.buf.partcs[ipart].rrad) < 0) {
      BMM_TLE_STDS();

      break;
    }

  if (fclose(stream) != 0) {
    BMM_TLE_STDS();

    return false;
  }

  return true;
}

static bool more_heresy(struct bmm_sdl const* const sdl) {
  FILE* const stream = fopen("more-heresy.data", "w");
  if (stream == NULL) {
    BMM_TLE_STDS();

    return false;
  }

  for (size_t ipart = 0; ipart < sdl->dem.buf.npart; ++ipart) {
    for (size_t ilink = 0; ilink < sdl->dem.buf.links[ipart].n; ++ilink) {
      size_t const jpart = sdl->dem.buf.links[ipart].linkl[ilink].i;

      if (fprintf(stream, "%zu %g %g\n",
            (size_t) 0,
            sdl->dem.buf.parts[ipart].lin.r[0],
            sdl->dem.buf.parts[ipart].lin.r[1]) < 0) {
        BMM_TLE_STDS();

        break; // out
      }

      if (fprintf(stream, "%zu %g %g\n",
            (size_t) 1,
            sdl->dem.buf.parts[jpart].lin.r[0],
            sdl->dem.buf.parts[jpart].lin.r[1]) < 0) {
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

static struct {
  double r;
  char name[4];
} const bmm_msg_covr[] = {
  {.r = 0.32, .name = "H"},
  {.r = 0.71, .name = "Ne"},
  {.r = 0.72, .name = "F"},
  {.r = 0.73, .name = "O"},
  {.r = 0.75, .name = "N"},
  {.r = 0.77, .name = "C"},
  {.r = 0.82, .name = "B"},
  {.r = 0.9, .name = "Be"},
  {.r = 0.93, .name = "He"},
  {.r = 0.98, .name = "Ar"},
  {.r = 0.99, .name = "Cl"},
  {.r = 1.02, .name = "S"},
  {.r = 1.06, .name = "P"},
  {.r = 1.11, .name = "Si"},
  {.r = 1.12, .name = "Kr"},
  {.r = 1.14, .name = "Br"},
  {.r = 1.15, .name = "Ni"},
  {.r = 1.16, .name = "Se"},
  {.r = 1.16, .name = "Co"},
  {.r = 1.17, .name = "Cu"},
  {.r = 1.17, .name = "Fe"},
  {.r = 1.17, .name = "Mn"},
  {.r = 1.18, .name = "Al"},
  {.r = 1.18, .name = "Cr"},
  {.r = 1.2, .name = "As"},
  {.r = 1.22, .name = "Ge"},
  {.r = 1.22, .name = "V"},
  {.r = 1.23, .name = "Li"},
  {.r = 1.25, .name = "Rh"},
  {.r = 1.25, .name = "Ru"},
  {.r = 1.25, .name = "Zn"},
  {.r = 1.26, .name = "Ga"},
  {.r = 1.26, .name = "Os"},
  {.r = 1.27, .name = "Ir"},
  {.r = 1.27, .name = "Tc"},
  {.r = 1.28, .name = "Re"},
  {.r = 1.28, .name = "Pd"},
  {.r = 1.3, .name = "W"},
  {.r = 1.3, .name = "Pt"},
  {.r = 1.3, .name = "Mo"},
  {.r = 1.31, .name = "Xe"},
  {.r = 1.32, .name = "Ti"},
  {.r = 1.33, .name = "I"},
  {.r = 1.34, .name = "Ta"},
  {.r = 1.34, .name = "Nb"},
  {.r = 1.34, .name = "Ag"},
  {.r = 1.34, .name = "Au"},
  {.r = 1.36, .name = "Te"},
  {.r = 1.36, .name = "Mg"},
  {.r = 1.41, .name = "Sn"},
  {.r = 1.41, .name = "Sb"},
  {.r = 1.42, .name = "U"},
  {.r = 1.44, .name = "In"},
  {.r = 1.44, .name = "Sc"},
  {.r = 1.44, .name = "Hf"},
  {.r = 1.45, .name = "Zr"},
  {.r = 1.45, .name = "At"},
  {.r = 1.46, .name = "Bi"},
  {.r = 1.46, .name = "Po"},
  {.r = 1.47, .name = "Pb"},
  {.r = 1.48, .name = "Cd"},
  {.r = 1.48, .name = "Tl"},
  {.r = 1.49, .name = "Hg"},
  {.r = 1.54, .name = "Na"},
  {.r = 1.56, .name = "Tm"},
  {.r = 1.56, .name = "Lu"},
  {.r = 1.57, .name = "Er"},
  {.r = 1.58, .name = "Ho"},
  {.r = 1.59, .name = "Dy"},
  {.r = 1.59, .name = "Tb"},
  {.r = 1.61, .name = "Gd"},
  {.r = 1.62, .name = "Y"},
  {.r = 1.62, .name = "Sm"},
  {.r = 1.63, .name = "Pm"},
  {.r = 1.64, .name = "Nd"},
  {.r = 1.65, .name = "Th"},
  {.r = 1.65, .name = "Ce"},
  {.r = 1.65, .name = "Pr"},
  {.r = 1.69, .name = "La"},
  {.r = 1.74, .name = "Yb"},
  {.r = 1.74, .name = "Ca"},
  {.r = 1.85, .name = "Eu"},
  {.r = 1.91, .name = "Sr"},
  {.r = 1.98, .name = "Ba"},
  {.r = 2.03, .name = "K"},
  {.r = 2.16, .name = "Rb"},
  {.r = 2.35, .name = "Cs"},
};

static bool serious_heresy(struct bmm_sdl const* const sdl) {
  FILE* const stream = fopen("heresy.xyz", "w");
  if (stream == NULL) {
    BMM_TLE_STDS();

    return false;
  }

  if (fprintf(stream, "%zu\n.\n", sdl->dem.buf.npart) < 0) {
    BMM_TLE_STDS();

    return false;
  }

  double const sfact = sdl->dem.opts.rmean;

  for (size_t ipart = 0; ipart < sdl->dem.buf.npart; ++ipart)
    if (fprintf(stream, "S %g %g 0.0 %g\n",
          sdl->dem.buf.parts[ipart].lin.r[0] / sfact,
          sdl->dem.buf.parts[ipart].lin.r[1] / sfact,
          sdl->dem.buf.partcs[ipart].rrad / sfact) < 0) {
      BMM_TLE_STDS();

      break;
    }

  if (fclose(stream) != 0) {
    BMM_TLE_STDS();

    return false;
  }

  return true;
}

static bool bmm_sdl_work(struct bmm_sdl* const sdl) {
  if (!bmm_sdl_video(sdl, sdl->width, sdl->height))
    return false;

  glClearColor4fv(glBlack);
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);

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
                  (double) sdl->width * 0.5, (double) sdl->height * 0.5, 0.5);
              break;
            case SDLK_PLUS:
            case SDLK_KP_PLUS:
              bmm_sdl_zoom(sdl,
                  (double) sdl->width * 0.5, (double) sdl->height * 0.5, 2.0);
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
              sdl->itarget = (size_t) rand() % sdl->dem.buf.npart;
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
  }
}

bool bmm_sdl_run(struct bmm_sdl* const sdl) {
  if (SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 8) == -1 ||
      SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE, 8) == -1 ||
      SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE, 8) == -1 ||
      SDL_GL_SetAttribute(SDL_GL_ALPHA_SIZE, 8) == -1 ||
      SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 8) == -1 ||
      SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1) == -1) {
    BMM_TLE_EXTS(BMM_TLE_SDL, "SDL error: %s", SDL_GetError());

    return false;
  }

  if (sdl->opts.ms > 0)
    if (SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS, 1) == -1 ||
        SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES,
          (int) sdl->opts.ms) == -1) {
      BMM_TLE_EXTS(BMM_TLE_SDL, "SDL error: %s", SDL_GetError());

      return false;
    }

  // TODO This is a bad idea and should not even work.
  int argc = 1;
  char arg[] = "";
  char* argv[] = {arg};
  glutInit(&argc, argv);

  if (!bmm_sdl_work(sdl))
    return false;

  return true;
}

static bool bmm_sdl_run_sdl(struct bmm_sdl_opts const* const opts) {
  struct bmm_sdl* const sdl = malloc(sizeof *sdl);
  if (sdl == NULL) {
    BMM_TLE_STDS();

    return false;
  }

  bmm_sdl_def(sdl, opts);
  bool const result = bmm_sdl_run(sdl);

  free(sdl);

  return result;
}

bool bmm_sdl_run_with(struct bmm_sdl_opts const* const opts) {
  if (SDL_Init(SDL_INIT_VIDEO) == -1) {
    BMM_TLE_EXTS(BMM_TLE_SDL, "SDL error: %s", SDL_GetError());

    return false;
  }

  bool const result = bmm_sdl_run_sdl(opts);

  SDL_Quit();

  return result;
}
