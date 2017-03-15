#include "fp.h"
#include "gl.h"
#include <GL/gl.h>
#include <math.h>
#include <stddef.h>

extern inline void glClearColor4fv(GLfloat const*);

extern inline void glClearColor3fv(GLfloat const*);

void glSkewedAnnulus(GLfloat const x, GLfloat const y,
    GLfloat const r, GLfloat const rhole, GLfloat const rskew,
    GLfloat const askew, size_t const ncorner, GLfloat const* const c4fv) {
  glBegin(GL_TRIANGLE_STRIP);

  glColor4fv(c4fv);

  GLfloat const xskew = x + rskew * cosf(askew);
  GLfloat const yskew = y + rskew * sinf(askew);
  GLfloat const acorner = (float) M_2PI / (float) ncorner;

  glVertex2f(xskew + rhole, yskew);
  glVertex2f(x + r, y);
  for (size_t i = 1; i < ncorner; ++i) {
    GLfloat const a = (float) i * acorner;
    GLfloat const cosa = cosf(a);
    GLfloat const sina = sinf(a);

    glVertex2f(xskew + rhole * cosa, yskew + rhole * sina);
    glVertex2f(x + r * cosa, y + r * sina);
  }
  glVertex2f(xskew + rhole, yskew);
  glVertex2f(x + r, y);

  glEnd();
}

void glAnnulus(GLfloat const x, GLfloat const y,
    GLfloat const r, GLfloat const rhole,
    size_t const ncorner, GLfloat const* const c4fv) {
  glBegin(GL_TRIANGLE_STRIP);

  glColor4fv(c4fv);

  GLfloat const acorner = (float) M_2PI / (float) ncorner;

  glVertex2f(x + rhole, y);
  glVertex2f(x + r, y);
  for (size_t i = 1; i < ncorner; ++i) {
    GLfloat const a = (float) i * acorner;
    GLfloat const cosa = cosf(a);
    GLfloat const sina = sinf(a);

    glVertex2f(x + rhole * cosa, y + rhole * sina);
    glVertex2f(x + r * cosa, y + r * sina);
  }
  glVertex2f(x + rhole, y);
  glVertex2f(x + r, y);

  glEnd();
}

void glDisc(GLfloat const x, GLfloat const y, GLfloat const r,
    size_t const ncorner, GLfloat const* const c4fv) {
  glBegin(GL_TRIANGLE_FAN);

  glColor4fv(c4fv);

  GLfloat const acorner = (float) M_2PI / (float) ncorner;

  glVertex2f(x, y);
  glVertex2f(x + r, y);
  for (size_t i = 1; i < ncorner; ++i) {
    GLfloat const a = (float) i * acorner;
    GLfloat const cosa = cosf(a);
    GLfloat const sina = sinf(a);

    glVertex2f(x + r * cosa, y + r * sina);
  }
  glVertex2f(x + r, y);

  glEnd();
}
