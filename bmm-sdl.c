#include "dem.h"
#include "err.h"
#include "ext.h"
#include "msg.h"
#include "io.h"
#include <GL/gl.h>
#include <SDL/SDL.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

// TODO All of this.

void glFilledCircle(GLfloat const x, GLfloat const y, GLfloat const r,
    size_t const n) {
  static GLfloat const twopi = 2.0f * M_PI;

  glBegin(GL_TRIANGLE_FAN);

  glVertex2f(x, y);
  glVertex2f(x + r, y);
  for (size_t i = 1; i < n; ++i)
    glVertex2f(x + r * cos(i * twopi / n), y + r * sin(i * twopi / n));
  glVertex2f(x + r, y);

  glEnd();
}

static void draw_screen(struct bmm_state const* const state) {
  static float angle = 0.0f;
  static GLfloat v0[] = { -1.0f, -1.0f, 0.0f };
  static GLfloat v1[] = { 1.0f, -1.0f, 0.0f };
  static GLfloat v2[] = { 0.0f, 0.7f, 0.0f };
  static GLfloat white[] = { 1.0f, 1.0f, 1.0f };
  static GLfloat red[] = { 1.0f, 0.0f, 0.0f };
  static GLfloat green[] = { 0.0f, 1.0f, 0.0f };
  static GLfloat blue[] = { 0.0f, 0.0f, 1.0f };

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glScalef(40.0f, 40.0f, 1.0f);

  glTranslatef(7.0f, 5.0f, 0.0f);

  glRotatef(angle, 0.0f, 0.0f, 1.0f);
  if (++angle > 360.0f)
    angle = 0.0f;

  glBegin(GL_TRIANGLES);

  glColor3fv(red);
  glVertex3fv(v0);
  glColor3fv(green);
  glVertex3fv(v1);
  glColor3fv(blue);
  glVertex3fv(v2);

  glEnd();

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glColor3fv(white);
  size_t ipart = 0; // for (size_t ipart = 0; ipart < PART_MAX; ++ipart)
  glFilledCircle((float) state->parts[ipart].rpos[0],
      (float) state->parts[ipart].rpos[1],
      (float) state->parts[ipart].rrad + 10.0f,
      16);

  SDL_GL_SwapBuffers();
}

__attribute__((__nonnull__))
int main(int const argc, char **const argv) {
  if (SDL_Init(SDL_INIT_TIMER | SDL_INIT_VIDEO) == -1) {
    BMM_ERR_FWARN(SDL_Init,
        "Failed to initialize SDL because '%s'.", SDL_GetError());

    return EXIT_FAILURE;
  }

  SDL_WM_SetCaption("BMM", NULL);

  SDL_VideoInfo const* info = SDL_GetVideoInfo();
  if (info == NULL) {
    BMM_ERR_FWARN(SDL_GetVideoInfo,
        "Failed to get video information because '%s'.", SDL_GetError());

    return EXIT_FAILURE;
  }

  if (SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 5) == -1 ||
      SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE, 5) == -1 ||
      SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE, 5) == -1 ||
      SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 16) == -1 ||
      SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1)) {
    BMM_ERR_FWARN(SDL_GL_SetAttribute,
        "Failed to set OpenGL attributes because '%s'.", SDL_GetError());

    return EXIT_FAILURE;
  }

  int const width = 640;
  int const height = 480;
  int const bpp = info->vfmt->BitsPerPixel;
  int const flags = SDL_OPENGL;
  if (SDL_SetVideoMode(width, height, bpp, flags) == NULL) {
    BMM_ERR_FWARN(SDL_SetVideoMode,
        "Failed to set video mode because '%s'.", SDL_GetError());

    return EXIT_FAILURE;
  }

  glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
  glViewport(0, 0, width, height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  // TODO Wrangle this.
  glOrtho(-(double) width, (double) width, (double) height, -(double) height, 1.0, -1.0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  // Monitor `stdin` and draw at `fps`.
  // This must work whether `stdin` is stalled or not.
  //
  // set initial state
  // mark system stalled
  // for ever
  //   if frame time has passed
  //     draw frame
  //     mark system stalled
  //     reset frame time
  //   else
  //     if system is stalled
  //       read message with timeout of remaining frame time
  //       if read succeeded and message is relevant
  //         update state
  //         unmark system stalled
  //     else
  //       sleep remaining frame time

  int const fps = 30;
  Uint32 const fticks = (Uint32) (1000 / fps);

  struct bmm_state state;
  bmm_defstate(&state);

  bool stalled = true;

  Uint32 ticks = 0;

  for ever {
    SDL_Event event;
    while (SDL_PollEvent(&event)) {
      switch (event.type) {
        case SDL_QUIT:
          goto hell;
        case SDL_KEYDOWN:
          switch (event.key.keysym.sym) {
            case SDLK_ESCAPE:
            case SDLK_q:
              goto hell;
          }
      }
    }

    Uint32 const sticks = SDL_GetTicks();
    Uint32 const dticks = sticks - ticks;

    if (dticks >= fticks) {
      draw_screen(&state);

      stalled = true;

      ticks += fticks;
    } else {
      if (stalled) {
        struct timeval timeout;
        timeout.tv_sec = 0;
        timeout.tv_usec = 1000000 / fps;

        struct bmm_head head;
        switch (bmm_io_waitin(&timeout)) {
          case BMM_IO_ERROR:
            goto hell;
          case BMM_IO_READY:
            bmm_msg_get(&head, &state);

            stalled = false;
        }
      } else
        SDL_Delay(dticks);
    }
  }

hell:
  SDL_Quit();

  return EXIT_SUCCESS;
}
