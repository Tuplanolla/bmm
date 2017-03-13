#include "err.h"
#include "ext.h"
#include "io.h"
#include <GL/gl.h>
#include <SDL/SDL.h>
#include <stdio.h>
#include <stdlib.h>

// TODO All of this.

static void draw_screen( void ) {
  static float angle = 0.0f;
  static GLfloat v0[] = { -1.0f, -1.0f, 0.0f };
  static GLfloat v1[] = { 1.0f, -1.0f, 0.0f };
  static GLfloat v2[] = { 0.0f, 0.7f, 0.0f };
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

  glColor4fv(red);
  glVertex3fv(v0);
  glColor4fv(green);
  glVertex3fv(v1);
  glColor4fv(blue);
  glVertex3fv(v2);

  glEnd();

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

  glOrtho(0.0, (double) width, (double) height, 0.0, 1.0, -1.0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  // Monitor `stdin` and draw at `fps`.
  // This must work whether `stdin` is stalled or not.
  //
  // set initial state
  // mark system stalled
  // for ever
  //   if frame time has passed
  //     preempt frame time
  //     draw frame
  //     mark system stalled
  //     reset preempted frame time
  //   else
  //     if system is stalled
  //       preempt monitor time
  //       read message with timeout of remaining frame time
  //       if read succeeded and message is relevant
  //         update state
  //         unmark system stalled
  //     else
  //       sleep remaining frame time

  Uint32 ticks = 0;
  for ever {
    ticks = SDL_GetTicks();
    draw_screen();

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

    // TODO This limits the frame rate below `fps`, not to `fps`.
    int const fps = 60;
    SDL_Delay(1000 / fps);
  }

hell:
  SDL_Quit();

  return EXIT_SUCCESS;
}
