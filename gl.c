#include "fp.h"
#include "gl.h"
#include <GL/gl.h>
#include <math.h>
#include <stddef.h>

void glDisc(GLfloat const x, GLfloat const y, GLfloat const r,
    size_t const n, GLfloat const* const c4fv) {
  glBegin(GL_TRIANGLE_FAN);

  glColor4fv(c4fv);
  glVertex2f(x, y);
  glVertex2f(x + r, y);
  for (size_t i = 1; i < n; ++i) {
    GLfloat const a = (float) M_2PI * (float) i / (float) n;
    glVertex2f(x + r * cosf(a), y + r * sinf(a));
  }
  glVertex2f(x + r, y);

  glEnd();
}

void glPointedDisc(GLfloat const x, GLfloat const y, GLfloat const r,
    GLfloat const a, size_t const n,
    GLfloat const* const c4fv, GLfloat const* const mc4fv) {
  glDisc(x, y, r, n, c4fv);

  GLfloat const r2 = 0.5f * r;
  GLfloat const r4 = 0.5f * r2;
  glDisc(x + r2 * cosf(a), y + r2 * sinf(a), r4, n, mc4fv);
}
