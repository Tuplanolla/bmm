#ifndef BMM_GL2_H
/// Utilities for OpenGL 2.0 and later.
#define BMM_GL2_H

#include <GL/glew.h>

#include <stdbool.h>

#include "ext.h"

__attribute__ ((__nonnull__))
void bmm_gl2_tle(GLuint, PFNGLGETSHADERIVPROC, PFNGLGETSHADERINFOLOGPROC);

__attribute__ ((__nonnull__))
bool bmm_gl2_compile(GLuint *, GLenum, char const *);

__attribute__ ((__nonnull__))
bool bmm_gl2_link(GLuint *, GLuint, GLuint);

__attribute__ ((__nonnull__))
bool bmm_gl2_build(GLuint *, GLuint *, GLuint *, char const *, char const *);

__attribute__ ((__nonnull__))
bool bmm_gl2_buffer(GLuint *, GLenum, GLenum, GLvoid const *, GLsizei);

#endif
