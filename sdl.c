#include "dem.h"
#include "err.h"
#include "ext.h"
#include "fp.h"
#include "gl.h"
#include "io.h"
#include "msg.h"
#include "sdl.h"
#include <GL/gl.h>
#include <SDL/SDL.h>
#include <stdbool.h>
#include <stddef.h>
#include <sys/time.h>

extern inline void bmm_sdl_t_to_timeval(struct timeval*, Uint32);

extern inline Uint32 bmm_sdl_t_from_timeval(struct timeval const*);

extern inline Uint32 bmm_sdl_trem(Uint32, Uint32);

void bmm_sdl_defopts(struct bmm_sdl_opts* const opts) {
  opts->width = 640;
  opts->height = 480;
  opts->fps = 20;
  opts->ms = 8;
}

void bmm_sdl_def(struct bmm_sdl* const sdl,
    struct bmm_sdl_opts const* const opts) {
  sdl->opts = *opts;
  sdl->width = opts->width;
  sdl->height = opts->height;
  sdl->qaspect = (double) opts->width / (double) opts->height;
  sdl->qzoom = 1.0;
  sdl->rorigin[0] = 0.0;
  sdl->rorigin[1] = 0.0;
  sdl->tstep = opts->fps > 1000 ? 1 : (Uint32) (1000 / opts->fps);
  sdl->pstale = true;

  struct bmm_dem_opts defopts;
  bmm_dem_defopts(&defopts);
  bmm_dem_def(&sdl->dem, &defopts);
}

// TODO All kinds of things.

static void bmm_sdl_proj(struct bmm_sdl const* const sdl,
    double* const xproj, double* const yproj,
    double* const wproj, double* const hproj) {
  double const w = sdl->dem.rexts[0];
  double const h = sdl->dem.rexts[1];
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
  *xproj = xorigin - (weither - wzoom) / 2.0;
  *yproj = yorigin - (heither - hzoom) / 2.0;
}

static void bmm_sdl_zoom(struct bmm_sdl* const sdl,
    int const xscreen, int const yscreen, double const q) {
  double xproj;
  double yproj;
  double wproj;
  double hproj;
  bmm_sdl_proj(sdl, &xproj, &yproj, &wproj, &hproj);

  double const x = bmm_fp_lerp((double) xscreen,
      0.0, (double) sdl->width, xproj, xproj + wproj);
  double const y = bmm_fp_lerp((double) yscreen,
      (double) sdl->height, 0.0, yproj, yproj + hproj);

  double const ox = sdl->rorigin[0];
  double const oy = sdl->rorigin[1];

  sdl->rorigin[0] = x;
  sdl->rorigin[1] = y;

  sdl->qzoom *= q;

  bmm_sdl_proj(sdl, &xproj, &yproj, &wproj, &hproj);

  double const x2 = bmm_fp_lerp((double) xscreen,
      0.0, (double) sdl->width, xproj, xproj + wproj);
  double const y2 = bmm_fp_lerp((double) yscreen,
      (double) sdl->height, 0.0, yproj, yproj + hproj);

  sdl->rorigin[0] -= x2 - x;
  sdl->rorigin[1] -= y2 - y;
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

  size_t const ncorner = 32;

  glColor4fv(sdl->pstale ? glRed : glGreen);
  glDisc(0.05f, 0.05f, 0.025f, ncorner);

  glColor4fv(glYellow);
  for (size_t ipart = 0; ipart < BMM_PART_MAX; ++ipart) {
    float const x = (float) sdl->dem.parts[ipart].rpos[0] + 0.5f;
    float const y = (float) sdl->dem.parts[ipart].rpos[1] + 0.5f;
    float const r = (float) sdl->dem.parts[ipart].rrad + 0.1f;

    glSkewedAnnulus(x, y, r, r / 4.0f, r / 2.0f, 0.0f, ncorner);
  }

  glColor4fv(glWhite);
  glRectWire(0.0f, 0.0f, (float) sdl->dem.rexts[0], (float) sdl->dem.rexts[1]);

  SDL_GL_SwapBuffers();
}

static bool bmm_sdl_work(struct bmm_sdl* const sdl) {
  Uint32 tnow = SDL_GetTicks();
  Uint32 tnext = tnow + sdl->tstep;
  Uint32 trem = bmm_sdl_trem(tnow, tnext);

  for ever {
    // Handle events just before drawing for maximal responsiveness.
    SDL_Event event;
    while (SDL_PollEvent(&event))
      switch (event.type) {
        case SDL_QUIT:
          return true;
        case SDL_VIDEORESIZE:
          { // TODO Say no to copy-paste.
          SDL_VideoInfo const* const info = SDL_GetVideoInfo();
          if (info == NULL) {
            BMM_ERR_FWARN(SDL_GetVideoInfo, "SDL error: %s", SDL_GetError());

            return false;
          }

          int const width = event.resize.w;
          int const height = event.resize.h;
          int const bpp = info->vfmt->BitsPerPixel;
          int const flags = SDL_OPENGL | SDL_RESIZABLE;
          if (SDL_SetVideoMode(width, height, bpp, flags) == NULL) {
            BMM_ERR_FWARN(SDL_SetVideoMode, "SDL error: %s", SDL_GetError());

            return false;
          }

          glClearColor4fv(glBlack);
          glViewport(0, 0, width, height);
          glEnable(GL_CULL_FACE);
          glCullFace(GL_BACK);

          sdl->width = (unsigned int) width;
          sdl->height = (unsigned int) height;
          sdl->qaspect = (double) width / (double) height;
          }
          break;
        case SDL_KEYDOWN:
          switch (event.key.keysym.sym) {
            case SDLK_ESCAPE:
            case SDLK_q:
              return true;
            case SDLK_PLUS:
            case SDLK_KP_PLUS:
              sdl->qzoom *= 2.0;
              break;
            case SDLK_MINUS:
            case SDLK_KP_MINUS:
              sdl->qzoom *= 0.5;
              break;
            case SDLK_LEFT:
              sdl->rorigin[0] -= 0.25;
              break;
            case SDLK_RIGHT:
              sdl->rorigin[0] += 0.25;
              break;
            case SDLK_DOWN:
              sdl->rorigin[1] -= 0.25;
              break;
            case SDLK_UP:
              sdl->rorigin[1] += 0.25;
              break;
            case SDLK_0:
            case SDLK_KP0:
              sdl->qzoom = 1.0;
              sdl->rorigin[0] = 0.0;
              sdl->rorigin[1] = 0.0;
              break;
          }
          break;
        case SDL_MOUSEBUTTONDOWN:
          switch (event.button.button) {
            case SDL_BUTTON_WHEELDOWN:
              bmm_sdl_zoom(sdl, event.button.x, event.button.y, 0.5);
              break;
            case SDL_BUTTON_WHEELUP:
              bmm_sdl_zoom(sdl, event.button.x, event.button.y, 2.0);
              break;
          }
      }

    // Draw before state transformations to account for the initial state.
    bmm_sdl_draw(sdl);

    // Recompute the remaining time after event handling and drawing
    // in case they take a while to complete.
    tnow = SDL_GetTicks();
    trem = bmm_sdl_trem(tnow, tnext);

    // Use the remaining time to wait for input.
    struct timeval timeout;
    bmm_sdl_t_to_timeval(&timeout, trem);
    struct bmm_head head;
    switch (bmm_io_waitin(&timeout)) {
      case BMM_IO_ERROR:
        return false;
      case BMM_IO_READY:
        if (bmm_msg_get(&head, &sdl->dem))
          sdl->pstale = false;
        else
      case BMM_IO_TIMEOUT:
          sdl->pstale = true;
    }

    // Recompute the remaining time again after waiting for input
    // in case it arrives early or takes a while to process.
    tnow = SDL_GetTicks();
    trem = bmm_sdl_trem(tnow, tnext);

    // Sleep the remaining time and allow the tick counter to wrap.
    SDL_Delay(trem);
    tnext += sdl->tstep;
  }
}

bool bmm_sdl_run(struct bmm_sdl_opts const* const opts) {
  bool result = true;

  if (SDL_Init(SDL_INIT_VIDEO) == -1) {
    BMM_ERR_FWARN(SDL_Init, "SDL error: %s", SDL_GetError());

    goto ret_error;
  }

  SDL_WM_SetCaption("BMM", NULL);

  if (SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 8) == -1 ||
      SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE, 8) == -1 ||
      SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE, 8) == -1 ||
      SDL_GL_SetAttribute(SDL_GL_ALPHA_SIZE, 8) == -1 ||
      SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 8) == -1 ||
      SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1) == -1) {
    BMM_ERR_FWARN(SDL_GL_SetAttribute, "SDL error: %s", SDL_GetError());

    goto quit_error;
  }

  if (opts->ms > 0)
    if (SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS, 1) == -1 ||
        SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, (int) opts->ms) == -1) {
      BMM_ERR_FWARN(SDL_GL_SetAttribute, "SDL error: %s", SDL_GetError());

      goto quit_error;
    }

  SDL_VideoInfo const* const info = SDL_GetVideoInfo();
  if (info == NULL) {
    BMM_ERR_FWARN(SDL_GetVideoInfo, "SDL error: %s", SDL_GetError());

    goto quit_error;
  }

  int const width = (int) opts->width;
  int const height = (int) opts->height;
  int const bpp = info->vfmt->BitsPerPixel;
  int const flags = SDL_OPENGL | SDL_RESIZABLE;
  if (SDL_SetVideoMode(width, height, bpp, flags) == NULL) {
    BMM_ERR_FWARN(SDL_SetVideoMode, "SDL error: %s", SDL_GetError());

    goto quit_error;
  }

  glClearColor4fv(glBlack);
  glViewport(0, 0, width, height);
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);

  struct bmm_sdl sdl;
  bmm_sdl_def(&sdl, opts);
  if (!bmm_sdl_work(&sdl))
    goto quit_error;

  goto quit;
quit_error:
  result = false;
quit:
  SDL_Quit();

  goto ret;
ret_error:
  result = false;
ret:
  return result;
}
