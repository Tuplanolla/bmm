#ifndef BMM_GLUT_H
/// Interactive viewer based on OpenGL and FreeGLUT.
#define BMM_GLUT_H

#include <GL/glew.h>

#include <setjmp.h>
#include <stdbool.h>
#include <stddef.h>

#include "ext.h"

/// This structure contains graphics options.
struct bmm_glut_opts {
  char const* vpath;
  char const* fpath;
};

/// This structure tracks resources.
struct bmm_glut {
  struct bmm_glut_opts opts;
  char* prog;
  jmp_buf env;
  GLuint vshader;
  GLuint fshader;
  GLuint program;
};

/// The call `bmm_glut_opts_def(opts)`
/// writes the default graphics options into `opts`.
/// All messages are stopped by default.
__attribute__ ((__nonnull__))
void bmm_glut_opts_def(struct bmm_glut_opts*);

/// The call `bmm_glut_def(glut)`
/// writes the default graphics state into `glut`.
__attribute__ ((__nonnull__))
void bmm_glut_def(struct bmm_glut*, struct bmm_glut_opts const*);

/// The call `bmm_glut_step(glut)`
/// processes one iglutoming message with the graphics state `glut`.
__attribute__ ((__nonnull__))
enum bmm_io_read bmm_glut_step(struct bmm_glut*);

/// The call `bmm_glut_run(glut)`
/// processes all iglutoming messages and
/// handles signals with the graphics state `glut`.
__attribute__ ((__nonnull__))
bool bmm_glut_run(struct bmm_glut*);

/// The call `bmm_glut_run_with(opts)`
/// processes all iglutoming messages and
/// handles signals with the graphics options `opts`.
__attribute__ ((__nonnull__))
bool bmm_glut_run_with(struct bmm_glut_opts const*);

#endif
