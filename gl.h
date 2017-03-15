// Graphics operations for OpenGL.
#ifndef BMM_GL_H
#define BMM_GL_H

#include <GL/gl.h>
#include <stddef.h>

// The call `glDisc(x, y, r, n, c4fv)`
// draws a filled circle centered at `x` and `y`
// with the radius `r`, segment number `n` and fill color `c4fv`.
void glDisc(GLfloat, GLfloat, GLfloat, size_t, GLfloat const*);

// The call `glPointedDisc(x, y, r, a, n, c4fv, mc4fv)`
// draws a rotatable filled circle centered at `x` and `y`
// with the radius `r`, angle `a`, segment number `n`,
// fill color `c4fv` and marker color `mc4fv`.
void glPointedDisc(GLfloat, GLfloat, GLfloat, GLfloat,
    size_t, GLfloat const*, GLfloat const*);

#endif
