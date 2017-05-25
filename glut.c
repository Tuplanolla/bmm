#include <GL/glew.h>
#ifdef FREE
#include <GL/freeglut.h>
#else
#include <GL/glut.h>
#endif

#include <math.h>
#include <setjmp.h>
#include <signal.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gl2.h"
#include "glut.h"

#include "conf.h"
// TODO Undepend.
#include "dem.h"
#include "ext.h"
#include "io.h"
#include "msg.h"
#include "sig.h"
#include "tle.h"

/// This is necessary since GLUT does not support passing in closures.
static struct bmm_glut* bmm_glut;

void bmm_glut_opts_def(struct bmm_glut_opts* const opts) {
  opts->vpath = "bmm-glut.vs.glsl";
  opts->fpath = "bmm-glut.fs.glsl";
}

void bmm_glut_def(struct bmm_glut* const glut,
    struct bmm_glut_opts const* const opts) {
  glut->opts = *opts;
  glut->prog = NULL;
  glut->vshader = 0;
  glut->fshader = 0;
  glut->program = 0;
}

static enum bmm_io_read msg_read(void* buf, size_t const n,
    __attribute__ ((__unused__)) void* const ptr) {
  return bmm_io_readin(buf, n);
}

static bool bmm_glut_open(struct bmm_glut* const glut) {

  return true;
}

static bool bmm_glut_close(struct bmm_glut* const glut) {

  return true;
}

enum bmm_io_read bmm_glut_step(struct bmm_glut* const glut) {
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

  enum bmm_msg_num num;
  switch (bmm_msg_num_read(&num, msg_read, NULL)) {
    case BMM_IO_READ_EOF:
      BMM_TLE_EXTS(BMM_TLE_NUM_IO, "Unexpected end");
    case BMM_IO_READ_ERROR:
      return BMM_IO_READ_ERROR;
  }

  switch (num) {
    case BMM_MSG_NUM_PARTS:
      {
        size_t npart;

        switch (msg_read(&npart, sizeof npart, NULL)) {
          case BMM_IO_READ_EOF:
            BMM_TLE_EXTS(BMM_TLE_NUM_IO, "Unexpected end");
          case BMM_IO_READ_ERROR:
            return BMM_IO_READ_ERROR;
        }

        double parts[BMM_MPART];

        switch (msg_read(parts, sizeof parts, NULL)) {
          case BMM_IO_READ_EOF:
            BMM_TLE_EXTS(BMM_TLE_NUM_IO, "Unexpected end");
          case BMM_IO_READ_ERROR:
            return BMM_IO_READ_ERROR;
        }
      }

      break;
    default:
      dynamic_assert(false, "Unsupported message number");
  }

  return BMM_IO_READ_SUCCESS;
}

static bool bmm_glut_run_(struct bmm_glut* const glut) {
  for ever {
    int signum;
    if (bmm_sig_use(&signum))
      switch (signum) {
        case SIGINT:
        case SIGQUIT:
        case SIGTERM:
        case SIGPIPE:
          BMM_TLE_EXTS(BMM_TLE_NUM_ASYNC, "Interrupted");

          return false;
      }

    switch (bmm_glut_step(glut)) {
      case BMM_IO_READ_ERROR:
        return false;
      case BMM_IO_READ_EOF:
        return true;
    }
  }
}

static void ifunc(void) {
}

static void render(void) {
}

static void reshape(int const w, int const h) {
}

static void bmm_glut_exit(void) {
  if (bmm_glut->program != 0)
    glDeleteProgram(bmm_glut->program);

  if (bmm_glut->fshader != 0)
    glDeleteShader(bmm_glut->fshader);

  if (bmm_glut->vshader != 0)
    glDeleteShader(bmm_glut->vshader);

  free(bmm_glut->prog);
}

static void kfunc(unsigned char const key, int const x, int const y) {
  switch (key) {
    case '\x1b':
    case 'q':
#ifdef FREE
      glutLeaveMainLoop();
#else
      bmm_glut_exit();

      longjmp(bmm_glut->env, 1);
#endif
  }
}

static void sfunc(int const key, int const x, int const y) {
}

static void mfunc(int const button, int const state, int const x, int const y) {
  // The button parameter is one of GLUT_LEFT_BUTTON, GLUT_MIDDLE_BUTTON, or GLUT_RIGHT_BUTTON.
  // The state parameter is either GLUT_UP or GLUT_DOWN.
}

#ifdef FREE
static void mwfunc(int const wheel, int const dir, int const x, int const y) {
  // There is the wheel number, direction is +/- 1, and x and y are the mouse coordinates.
}
#endif

bool bmm_glut_run(struct bmm_glut* const glut) {
  int const sigs[] = {SIGINT, SIGQUIT, SIGTERM, SIGPIPE};
  if (bmm_sig_register(sigs, nmembof(sigs)) != SIZE_MAX) {
    BMM_TLE_STDS();

    return false;
  }

  bmm_glut = glut;

  char const* const prog = bmm_tle_prog();
  size_t const size = strlen(prog) + 1;

  bmm_glut->prog = malloc(size);
  if (bmm_glut->prog == NULL) {
    BMM_TLE_STDS();

    return false;
  }

  (void) memcpy(bmm_glut->prog, prog, size);

  int argc = 1;
  char* argv[] = {bmm_glut->prog, NULL};
  glutInit(&argc, argv);

#ifdef FREE
  glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_GLUTMAINLOOP_RETURNS);
#endif

  glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
  glutInitWindowSize(640, 480);
  int const window = glutCreateWindow(BMM_SHORT);
  glutSetWindow(window);

  glutDisplayFunc(render);
  glutIdleFunc(ifunc);
  glutReshapeFunc(reshape);

  glewInit();
  if (!GLEW_VERSION_3_3) {
    BMM_TLE_EXTS(BMM_TLE_NUM_GL, "Failed to initialize OpenGL 3.3");

    return false;
  }

  glutKeyboardFunc(kfunc);
  glutSpecialFunc(sfunc);
  glutMouseFunc(mfunc);
#ifdef FREE
  glutMouseWheelFunc(mwfunc);
#endif

  if (setjmp(bmm_glut->env) == 0)
    glutMainLoop();

#ifdef FREE
  bmm_glut_exit();
#endif

  return true;
}

bool bmm_glut_run_with(struct bmm_glut_opts const* const opts) {
  struct bmm_glut glut;
  bmm_glut_def(&glut, opts);

  return bmm_glut_run(&glut);
}
