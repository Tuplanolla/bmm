#include <GL/glew.h>

#ifdef FREE
#include <GL/freeglut.h>
#else
#include <GL/glut.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "io.h"
#include "tle.h"

// Thank you, `www.opengl-tutorial.org`.

// Dirty.
static void show_info_log(GLuint object,
    PFNGLGETSHADERIVPROC glGet__iv,
    PFNGLGETSHADERINFOLOGPROC glGet__InfoLog) {
  GLint log_length;
  char *log;

  glGet__iv(object, GL_INFO_LOG_LENGTH, &log_length);
  log = malloc(log_length);
  glGet__InfoLog(object, log_length, NULL, log);
  fprintf(stderr, "%s", log);
  free(log);
}

static GLuint bmm_gl2_compile(GLenum const num, char const* path) {
  GLuint const shader = glCreateShader(num);
  if (shader == 0)
    return 0;

  size_t size;
  void* buf = bmm_io_conts(&size, path);
  if (buf == NULL)
    return 0;

  GLchar const* const str = buf;
  GLint const length = (GLint) size;
  glShaderSource(shader, 1, &str, &length);

  free(buf);

  glCompileShader(shader);

  // Dirty.
  GLint params;
  glGetShaderiv(shader, GL_COMPILE_STATUS, &params);
  if (params == 0) {
    fprintf(stderr, "Failed to compile %s:\n", path);
    show_info_log(shader, glGetShaderiv, glGetShaderInfoLog);
    glDeleteShader(shader);

    return 0;
  }

  return shader;
}

static GLuint bmm_gl2_link(GLuint const vshader, GLuint const fshader) {
  GLuint program = glCreateProgram();
  glAttachShader(program, vshader);
  glAttachShader(program, fshader);
  glLinkProgram(program);

  // Dirty.
  GLint params;
  glGetProgramiv(program, GL_LINK_STATUS, &params);
  if (params == 0) {
    fprintf(stderr, "Failed to link shader program:\n");
    show_info_log(program, glGetProgramiv, glGetProgramInfoLog);
    glDeleteProgram(program);

    return 0;
  }

  return program;
}

// Dirty.
static struct {
  GLuint vertex_buffer, element_buffer;
  GLuint textures[2];

  GLuint vertex_shader, fragment_shader, program;

  struct {
    GLint fade_factor;
  } uniforms;

  struct {
    GLint position;
  } attributes;

  GLfloat fade_factor;
} g_resources;

// Dirty.
static GLuint make_buffer(
    GLenum target,
    const void *buffer_data,
    GLsizei buffer_size) {
  GLuint buffer;
  glGenBuffers(1, &buffer);
  glBindBuffer(target, buffer);
  glBufferData(target, buffer_size, buffer_data, GL_STATIC_DRAW);

  return buffer;
}

// Dirty.
static const GLfloat g_vertex_buffer_data[] = {
  -1.0f, -1.0f,
  1.0f, -1.0f,
  -1.0f,  1.0f,
  1.0f,  1.0f
};
static const GLushort g_element_buffer_data[] = { 0, 1, 2, 3 };

// Dirty.
static int make_resources(void) {
  g_resources.vertex_buffer = make_buffer(
      GL_ARRAY_BUFFER,
      g_vertex_buffer_data,
      sizeof g_vertex_buffer_data);
  g_resources.element_buffer = make_buffer(
      GL_ELEMENT_ARRAY_BUFFER,
      g_element_buffer_data,
      sizeof g_element_buffer_data);

  g_resources.vertex_shader = bmm_gl2_compile(GL_VERTEX_SHADER, "glut.vs.glsl");
  if (g_resources.vertex_shader == 0)
    return 0;

  g_resources.fragment_shader = bmm_gl2_compile(GL_FRAGMENT_SHADER, "glut.fs.glsl");
  if (g_resources.fragment_shader == 0)
    return 0;

  g_resources.program = bmm_gl2_link(
      g_resources.vertex_shader,
      g_resources.fragment_shader);
  if (g_resources.program == 0)
    return 0;

  g_resources.uniforms.fade_factor
    = glGetUniformLocation(g_resources.program, "fade_factor");

  g_resources.attributes.position
    = glGetAttribLocation(g_resources.program, "position");

  return 1;
}

// Dirty.
static void render(void) {
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT);

  glUseProgram(g_resources.program);

  glUniform1f(g_resources.uniforms.fade_factor, g_resources.fade_factor);

  glBindBuffer(GL_ARRAY_BUFFER, g_resources.vertex_buffer);
  glVertexAttribPointer(
      g_resources.attributes.position,  /* attribute */
      2,                                /* size */
      GL_FLOAT,                         /* num */
      GL_FALSE,                         /* normalized? */
      sizeof(GLfloat)*2,                /* stride */
      NULL                          /* array buffer offset */
      );
  glEnableVertexAttribArray(g_resources.attributes.position);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g_resources.element_buffer);
  glDrawElements(
      GL_TRIANGLE_STRIP,  /* mode */
      4,                  /* count */
      GL_UNSIGNED_SHORT,  /* num */
      NULL            /* element array buffer offset */
      );

  glDisableVertexAttribArray(g_resources.attributes.position);

  glutSwapBuffers();
}

// Dirty.
static void update_fade_factor(void) {
  int milliseconds = glutGet(GLUT_ELAPSED_TIME);
  g_resources.fade_factor = sinf((float)milliseconds * 0.001f) * 0.5f + 0.5f;
  glutPostRedisplay();
}

static void reshape(int const w, int const h) {
}

static void before_exit(void) {
  // TODO Clean up here.
  puts("HUA!");
}

static void kfunc(unsigned char const key, int const x, int const y) {
  switch (key) {
    case '\x1b':
    case 'q':
#ifdef FREE
      glutLeaveMainLoop();
#else
      before_exit();

      exit(EXIT_SUCCESS);
#endif
  }
}

static void sfunc(int const key, int const x, int const y) {
}

static void mfunc(int const button, int const state, int const x, int const y) {
  // The button parameter is one of GLUT_LEFT_BUTTON, GLUT_MIDDLE_BUTTON, or GLUT_RIGHT_BUTTON.
  // The state parameter is either GLUT_UP or GLUT_DOWN.
}

static void mwfunc(int const wheel, int const dir, int const x, int const y) {
  // There is the wheel number, direction is +/- 1, and x and y are the mouse coordinates.
}

// Dirty.
int main(int argc, char** argv) {
  glutInit(&argc, argv);

#ifdef FREE
  glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_GLUTMAINLOOP_RETURNS);
#endif

  glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
  glutInitWindowSize(640, 480);
  glutCreateWindow("GLUT");

  glutDisplayFunc(render);
  glutIdleFunc(update_fade_factor);
  glutReshapeFunc(reshape);

  glewInit();
  // glDrawArraysInstanced needs 3.1
  if (!GLEW_VERSION_3_3) {
    fprintf(stderr, "OpenGL not available\n");

    return 1;
  }

  if (!make_resources()) {
    fprintf(stderr, "Failed to load resources\n");

    return 1;
  }

  glutKeyboardFunc(kfunc);
  glutSpecialFunc(sfunc);
  glutMouseFunc(mfunc);
#ifdef FREE
  glutMouseWheelFunc(mwfunc);
#endif

  glutMainLoop();

#ifdef FREE
  before_exit();

  return EXIT_SUCCESS;
#else
  return EXIT_FAILURE;
#endif
}
