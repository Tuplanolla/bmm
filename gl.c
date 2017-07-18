#include <GL/gl.h>
#include <math.h>
#include <stddef.h>

#include "fp.h"
#include "gl.h"

extern inline void glClearColor4fv(GLfloat const *);

extern inline void glClearColor3fv(GLfloat const *);

void glSkewedAnnulus(GLfloat const x, GLfloat const y,
    GLfloat const r, GLfloat const rhole, GLfloat const rskew,
    GLfloat const askew, size_t const ncorner) {
  glBegin(GL_TRIANGLE_STRIP);

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
    GLfloat const r, GLfloat const rhole, size_t const ncorner) {
  glBegin(GL_TRIANGLE_STRIP);

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

void glDisk(GLfloat const x, GLfloat const y,
    GLfloat const r, size_t const ncorner) {
  glBegin(GL_TRIANGLE_FAN);

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

void glRectWire(GLfloat const x, GLfloat const y,
    GLfloat const w, GLfloat const h) {
  glBegin(GL_LINE_LOOP);

  glVertex2d(x, y);
  glVertex2d(x + w, y);
  glVertex2d(x + w, y + h);
  glVertex2d(x, y + h);

  glEnd();
}
