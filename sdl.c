#include "dem.h"
#include "err.h"
#include "ext.h"
#include "fp.h"
#include "gl.h"
#include "msg.h"
#include "io.h"
#include "opt.h"
#include "sdl.h"
#include <GL/gl.h>
#include <SDL/SDL.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

extern inline void bmm_sdl_t_to_timeval(struct timeval*, Uint32);

extern inline Uint32 bmm_sdl_t_from_timeval(struct timeval const*);

extern inline Uint32 bmm_sdl_trem(Uint32, Uint32);

void bmm_sdl_defopts(struct bmm_sdl_opts* const opts) {
  opts->width = 640;
  opts->height = 480;
  opts->fps = 20;
  opts->ms = 8;
}

void bmm_sdl_defstate(struct bmm_sdl_state* const state,
    struct bmm_sdl_opts const* const opts) {
  struct bmm_dem_opts defopts;
  bmm_defopts(&defopts);
  bmm_dem_defstate(&state->state, &defopts);

  state->tstep = opts->fps > 1000 ? 1 : (Uint32) (1000 / opts->fps);
  state->stale = true;
}

// TODO Choose the graphics coordinate system.

static void draw(struct bmm_sdl_state const* const state) {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  size_t const ncorner = 32;

  glDisc(-620.0f, -460.0f, 10.0f, ncorner, state->stale ? glRed : glGreen);

  for (size_t ipart = 0; ipart < BMM_PART_MAX; ++ipart) {
    float const x = (float) state->state.parts[ipart].rpos[0];
    float const y = (float) state->state.parts[ipart].rpos[1];
    float const r = (float) state->state.parts[ipart].rrad + 40.0f;

    glSkewedAnnulus(x, y, r, r / 4.0f, r / 2.0f, 0.0f, ncorner, glWhite);
  }

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

  { // TODO Wrangle OpenGL options.
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    glFrontFace(GL_CW);

    glClearColor4fv(glBlack);
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    glOrtho(-(double) width, (double) width,
        (double) height, -(double) height, 1.0, -1.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
  }

  struct bmm_sdl_state state;
  bmm_sdl_defstate(&state, opts);

  // Start timing.
  Uint32 tnow = SDL_GetTicks();
  Uint32 tnext = tnow + state.tstep;
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
          }
      }

    // Draw before state transformations to account for the initial state.
    draw(&state);

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
        if (bmm_msg_get(&head, &state.state))
          state.stale = false;
        else
      case BMM_IO_TIMEOUT:
          state.stale = true;
    }

    // Recompute the remaining time again after waiting for input
    // in case it arrives early or takes a while to process.
    tnow = SDL_GetTicks();
    trem = bmm_sdl_trem(tnow, tnext);

    // Sleep and allow the tick counter to wrap.
    SDL_Delay(trem);
    tnext += state.tstep;
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
