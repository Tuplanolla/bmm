#include <GL/gl.h>
#include <GL/glut.h>
#include <SDL2/SDL.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/time.h>

#include "dem.h"
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

void bmm_sdl_opts_def(struct bmm_sdl_opts* const opts) {
  opts->width = 640;
  opts->height = 480;
  opts->fps = 20;
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

  struct bmm_dem_buf const* const buf = bmm_dem_getrbuf(&sdl->dem);

  // Particles.
  for (size_t ipart = 0; ipart < buf->npart; ++ipart) {
    float const x = (float) buf->parts[ipart].lin.r[0];
    float const y = (float) buf->parts[ipart].lin.r[1];
    float const r = (float) buf->partcs[ipart].rrad;
    float const a = (float) buf->parts[ipart].ang.alpha;

    glColor4fv(buf->partcs[ipart].free ? glYellow : glWhite);

    glSkewedAnnulus(x, y, r, r * 0.25f, r * 0.5f, a, ncorner);
  }

  // Focus.
  if (sdl->itarget != SIZE_MAX) {
    glColor4fv(glBlue);
    size_t const ipart = sdl->itarget;

    // Neighbor markers.
    glBegin(GL_LINES);
    for (size_t ineigh = 0; ineigh < bmm_dem_sizey(&buf->neigh.neighs[ipart]); ++ineigh) {
      size_t const jpart = bmm_dem_gety(&buf->neigh.neighs[ipart], ineigh);

      bool p = true;
      for (size_t idim = 0; idim < 2; ++idim)
        p = p &&
          fabs(buf->parts[ipart].lin.r[idim] - buf->parts[jpart].lin.r[idim]) <
          sdl->dem.rext[idim] / 2.0;

      if (p) {
        glVertex2dv(buf->parts[ipart].lin.r);
        glVertex2dv(buf->parts[jpart].lin.r);
      }
    }
    glEnd();
  }

  // Structural links.
  glColor4fv(glMagenta);
  for (size_t ipart = 0; ipart < buf->npart; ++ipart) {
    glBegin(GL_LINES);
    for (size_t ineigh = 0; ineigh < bmm_dem_sizel(&buf->links[ipart]); ++ineigh) {
      size_t const jpart = bmm_dem_getl(&buf->links[ipart], ineigh);

      bool p = true;
      for (size_t idim = 0; idim < 2; ++idim)
        p = p &&
          fabs(buf->parts[ipart].lin.r[idim] - buf->parts[jpart].lin.r[idim]) <
          sdl->dem.rext[idim] / 2.0;

      if (p) {
        glVertex2dv(buf->parts[ipart].lin.r);
        glVertex2dv(buf->parts[jpart].lin.r);
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
      buf->neigh.tnext);
  glString(strbuf, 8, 8 + 15 * ioff++, glWhite, GLUT_BITMAP_9_BY_15);
  (void) snprintf(strbuf, sizeof strbuf, "e (coeff. of restit.) = %g",
      bmm_dem_cor(&sdl->dem));
  glString(strbuf, 8, 8 + 15 * ioff++, glWhite, GLUT_BITMAP_9_BY_15);
  (void) snprintf(strbuf, sizeof strbuf, "n (number of particles) = %zu",
      buf->npart);
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
              return true;
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
              {
                struct bmm_dem_buf const* const buf = bmm_dem_getrbuf(&sdl->dem);
                sdl->itarget = (size_t) rand() % buf->npart;
              }
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
    // TODO Currently this sleeps too much.

    if (sdl->active) {
      // Use the remaining time to wait for input.
      // If there is no time left,
      // try to read just one message to prevent congestion.
      struct timeval timeout;
again:
      bmm_sdl_t_to_timeval(&timeout, trem);
      struct bmm_msg_head head;
      switch (bmm_io_waitin(&timeout)) {
        case BMM_IO_WAIT_ERROR:
          return false;
        case BMM_IO_WAIT_READY:
          if (bmm_msg_get(&head, &sdl->dem))
            sdl->stale = false;
          else {
            trem = bmm_sdl_t_from_timeval(&timeout);
            if (trem > 0)
              goto again;
            else
              sdl->stale = true;
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
    // TODO Fix the effect of `filter` on the vfps.
    // tnext += bmm_sdl_tstep(sdl);
    tnext += (bmm_sdl_tstep(sdl) - 1) / 2 + 1;
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
