#ifndef BMM_GL_H
/// Additional operations for OpenGL.
#define BMM_GL_H

#include <GL/gl.h>
#include <stddef.h>

#include "ext.h"

/// Primaries based on ITU-R BT.601 weights.
static GLfloat const glBlack[] = {0.0f, 0.0f, 0.0f};
static GLfloat const glRed[] = {0.701f, 0.0f, 0.0f};
static GLfloat const glGreen[] = {0.0f, 0.413f, 0.0f};
static GLfloat const glBlue[] = {0.0f, 0.0f, 0.886f};
static GLfloat const glCyan[] = {0.0f, 0.6495f, 0.6495f};
static GLfloat const glMagenta[] = {0.7935f, 0.0f, 0.7935f};
static GLfloat const glYellow[] = {0.557f, 0.557f, 0.0f};
static GLfloat const glWhite[] = {0.6666667f, 0.6666667f, 0.6666667f};

/// The call `glClearColor4fv(v)` is equivalent
/// to `glClearColor(v[0], v[1], v[2], v[3])`.
__attribute__ ((__nonnull__))
inline void glClearColor4fv(GLfloat const* const v) {
  glClearColor(v[0], v[1], v[2], v[3]);
}

/// The call `glClearColor3fv(v)` is equivalent
/// to `glClearColor(v[0], v[1], v[2], 1.0f)`.
__attribute__ ((__nonnull__))
inline void glClearColor3fv(GLfloat const* const v) {
  glClearColor(v[0], v[1], v[2], 1.0f);
}

/// The call `glSkewedAnnulus(x, y, r, rhole, rskew, askew, ncorner, c4fv)`
/// draws a filled circle with a hole in it.
/// The circle itself has a radius `r` and is centered at `x` and `y`
/// while the hole has a radius of `rhole` and
/// is at an angle of `a` and distance of `rskew` from the center.
/// The shape is drawn with `ncorner` corners.
void glSkewedAnnulus(GLfloat, GLfloat, GLfloat, GLfloat, GLfloat, GLfloat,
    size_t);

/// The call `glAnnulus(x, y, r, rhole, ncorner)` is equivalent
/// to `glSkewedAnnulus(x, y, r, rhole, 0.0f, 0.0f, ncorner)`.
void glAnnulus(GLfloat, GLfloat, GLfloat, GLfloat, size_t);

/// The call `glDisk(x, y, r, ncorner)` is equivalent
/// to `glAnnulus(x, y, r, 0.0f, ncorner)`.
void glDisk(GLfloat, GLfloat, GLfloat, size_t);

/// The call `glRectWire(x, y, w, h)`
/// draws the outline of a rectangle
/// with a corner at `x` and `y`, a width of `w` and a height of `h`.
void glRectWire(GLfloat, GLfloat, GLfloat, GLfloat);

#endif
