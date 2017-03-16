#include "dem.h"
#include "err.h"
#include "ext.h"
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
  sdl->ratio = (double) opts->width / (double) opts->height;
  sdl->zoom = 1.0;
  sdl->focus[0] = 0.5;
  sdl->focus[1] = 0.5;
  sdl->tstep = opts->fps > 1000 ? 1 : (Uint32) (1000 / opts->fps);
  sdl->stale = true;

  struct bmm_dem_opts defopts;
  bmm_dem_defopts(&defopts);
  bmm_dem_def(&sdl->dem, &defopts);
}

// TODO Check this math.

static void draw(struct bmm_sdl const* const sdl) {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

  double const wexts = sdl->dem.rexts[0] / sdl->zoom;
  double const hexts = sdl->dem.rexts[1] / sdl->zoom;
  double const wrat = hexts * sdl->ratio;

  bool const p = wrat > wexts;

  double const w = p ? wrat : wexts;
  double const h = p ? hexts : wexts / sdl->ratio;

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(sdl->focus[0] - w / 2.0, sdl->focus[0] + w / 2.0,
      sdl->focus[1] - h / 2.0, sdl->focus[1] + h / 2.0,
      -1.0, 1.0);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glColor4fv(glBlue);
  glRectd(0.0, 0.0, sdl->dem.rexts[0], sdl->dem.rexts[1]);

  size_t const ncorner = 32;

  glColor4fv(sdl->stale ? glRed : glGreen);
  glDisc(0.05f, 0.05f, 0.025f, ncorner);

  glColor4fv(glYellow);
  for (size_t ipart = 0; ipart < BMM_PART_MAX; ++ipart) {
    float const x = (float) sdl->dem.parts[ipart].rpos[0] + 0.5f;
    float const y = (float) sdl->dem.parts[ipart].rpos[1] + 0.5f;
    float const r = (float) sdl->dem.parts[ipart].rrad + 0.1f;

    glSkewedAnnulus(x, y, r, r / 4.0f, r / 2.0f, 0.0f, ncorner);
  }

  glColor4fv(glWhite);
  glDisc((float) sdl->focus[0], (float) sdl->focus[1], 0.01f, ncorner);

  SDL_GL_SwapBuffers();
}

bool bmm_sdl_run(struct bmm_sdl_opts const* const opts) {
  bool result = true;

  if (SDL_Init(SDL_INIT_TIMER | SDL_INIT_VIDEO) == -1) {
    BMM_ERR_FWARN(SDL_Init, "SDL error: %s", SDL_GetError());

    goto ret_error;
  }

  SDL_WM_SetCaption("BMM", NULL);

  SDL_VideoInfo const* const info = SDL_GetVideoInfo();
  if (info == NULL) {
    BMM_ERR_FWARN(SDL_GetVideoInfo, "SDL error: %s", SDL_GetError());

    goto quit_error;
  }

  if (SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 5) == -1 ||
      SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE, 5) == -1 ||
      SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE, 5) == -1 ||
      SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 16) == -1 ||
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

  int const width = (int) opts->width;
  int const height = (int) opts->height;
  int const bpp = info->vfmt->BitsPerPixel;
  int const flags = SDL_OPENGL;
  if (SDL_SetVideoMode(width, height, bpp, flags) == NULL) {
    BMM_ERR_FWARN(SDL_SetVideoMode, "SDL error: %s", SDL_GetError());

    goto quit_error;
  }

  glClearColor4fv(glBlack);

  glViewport(0, 0, width, height);

  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);

  // Create the initial state once SDL and OpenGL have been initialized.
  struct bmm_sdl sdl;
  bmm_sdl_def(&sdl, opts);

  // Start timing.
  Uint32 tnow = SDL_GetTicks();
  Uint32 tnext = tnow + sdl.tstep;
  Uint32 trem = bmm_sdl_trem(tnow, tnext);

  for ever {
    // Handle events just before drawing for maximal responsiveness.
    SDL_Event event;
    while (SDL_PollEvent(&event))
      switch (event.type) {
        case SDL_QUIT:
          goto quit;
        case SDL_KEYDOWN:
          switch (event.key.keysym.sym) {
            case SDLK_ESCAPE:
            case SDLK_q:
              goto quit;
            case SDLK_PLUS:
            case SDLK_KP_PLUS:
              sdl.zoom *= 2.0;
              break;
            case SDLK_MINUS:
            case SDLK_KP_MINUS:
              sdl.zoom *= 0.5;
              break;
            case SDLK_LEFT:
              sdl.focus[0] -= 0.25;
              break;
            case SDLK_RIGHT:
              sdl.focus[0] += 0.25;
              break;
            case SDLK_DOWN:
              sdl.focus[1] -= 0.25;
              break;
            case SDLK_UP:
              sdl.focus[1] += 0.25;
              break;
          }
      }

    // Draw before state transformations to account for the initial state.
    draw(&sdl);

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
        goto quit_error;
      case BMM_IO_READY:
        if (bmm_msg_get(&head, &sdl.dem))
          sdl.stale = false;
        else
      case BMM_IO_TIMEOUT:
          sdl.stale = true;
    }

    // Recompute the remaining time again after waiting for input
    // in case it arrives early or takes a while to process.
    tnow = SDL_GetTicks();
    trem = bmm_sdl_trem(tnow, tnext);

    // Sleep and allow the tick counter to wrap.
    SDL_Delay(trem);
    tnext += sdl.tstep;
  }

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
