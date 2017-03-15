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
#include <stdio.h>
#include <stdlib.h>

extern inline void bmm_sdl_ticks_to_timeval(struct timeval*, Uint32);

extern inline Uint32 bmm_sdl_ticks_from_timeval(struct timeval const*);

void bmm_sdl_defopts(struct bmm_sdl_opts* const opts) {
  opts->width = 640;
  opts->height = 480;
  opts->fps = 20;
  opts->msaa = 8;
}

// TODO All of this.

static void draw_screen(struct bmm_state const* const state,
    bool const stalled) {
  static GLfloat black[] = {0.0f, 0.0f, 0.0f, 1.0f};
  static GLfloat white[] = {1.0f, 1.0f, 1.0f, 1.0f};
  static GLfloat red[] = {1.0f, 0.0f, 0.0f, 1.0f};
  static GLfloat green[] = {0.0f, 1.0f, 0.0f, 1.0f};
  static GLfloat blue[] = {0.0f, 0.0f, 1.0f, 1.0f};

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glDisc(-620.0f, -460.0f, 10.0f, 16, stalled ? red : green);

  for (size_t ipart = 0; ipart < BMM_PART_MAX; ++ipart)
    glPointedDisc((float) state->parts[ipart].rpos[0],
        (float) state->parts[ipart].rpos[1],
        (float) state->parts[ipart].rrad + 40.0f,
        0.0f,
        16, white, black);

  SDL_GL_SwapBuffers();
}

bool bmm_sdl_run(struct bmm_sdl_opts const* const opts) {
  if (SDL_Init(SDL_INIT_TIMER | SDL_INIT_VIDEO) == -1) {
    BMM_ERR_FWARN(SDL_Init, "SDL error: %s", SDL_GetError());

    goto ret;
  }

  SDL_WM_SetCaption("BMM", NULL);

  SDL_VideoInfo const* const info = SDL_GetVideoInfo();
  if (info == NULL) {
    BMM_ERR_FWARN(SDL_GetVideoInfo, "SDL error: %s", SDL_GetError());

    goto quit;
  }

  if (SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 5) == -1 ||
      SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE, 5) == -1 ||
      SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE, 5) == -1 ||
      SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 16) == -1 ||
      SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1) == -1) {
    BMM_ERR_FWARN(SDL_GL_SetAttribute, "SDL error: %s", SDL_GetError());

    goto quit;
  }

  if (opts->msaa > 0)
    if (SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS, 1) == -1 ||
        SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES,
          (int) opts->msaa) == -1) {
      BMM_ERR_FWARN(SDL_GL_SetAttribute, "SDL error: %s", SDL_GetError());

      goto quit;
    }

  int const width = (int) opts->width;
  int const height = (int) opts->height;
  int const bpp = info->vfmt->BitsPerPixel;
  int const flags = SDL_OPENGL;
  if (SDL_SetVideoMode(width, height, bpp, flags) == NULL) {
    BMM_ERR_FWARN(SDL_SetVideoMode, "SDL error: %s", SDL_GetError());

    goto quit;
  }

  glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
  glViewport(0, 0, width, height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  // TODO Wrangle this.
  glOrtho(-(double) width, (double) width,
      (double) height, -(double) height, 1.0, -1.0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  // Monitor `stdin` and draw at `fps`.
  // This must work whether `stdin` is stalled or not.

  Uint32 const tstep = (Uint32) (1000 / opts->fps);

  struct bmm_state state;
  bmm_defstate(&state);

  bool stalled = true;

  Uint32 tnext = SDL_GetTicks();

  for ever {
    SDL_Event event;
    while (SDL_PollEvent(&event)) {
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
    }

    Uint32 tnow = SDL_GetTicks();
    Uint32 trem = tnow < tnext ? tnext - tnow : 0;

    if (stalled) {
      struct timeval timeout;
      bmm_sdl_ticks_to_timeval(&timeout, trem);

      struct bmm_head head;
      switch (bmm_io_waitin(&timeout)) {
        case BMM_IO_ERROR:
          goto quit;
        case BMM_IO_READY:
          if (bmm_msg_get(&head, &state))
            stalled = false;
      }
    }

    tnow = SDL_GetTicks();
    trem = tnow < tnext ? tnext - tnow : 0;

    SDL_Delay(trem);

    draw_screen(&state, stalled);

    stalled = true;

    tnext += tstep;
  }

  bool p = true;

  goto quit_fine;
quit:
  p = false;
quit_fine:
  SDL_Quit();

  goto ret_fine;
ret:
  p = false;
ret_fine:
  return p;
}
