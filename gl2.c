#include <GL/glew.h>

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>

#include "gl2.h"

#include "io.h"
#include "tle.h"

void bmm_gl2_tle(GLuint const object,
    PFNGLGETSHADERIVPROC const iv,
    PFNGLGETSHADERINFOLOGPROC const infoLog) {
  GLint length;
  iv(object, GL_INFO_LOG_LENGTH, &length);
  if (length <= 0)
    return;

  size_t const size = (size_t) length;

  char *const buf = malloc(size);
  if (buf == NULL) {
    BMM_TLE_STDS();

    return;
  }

  infoLog(object, (GLsizei) size, NULL, (GLchar *) buf);
  BMM_TLE_EXTS(BMM_TLE_NUM_GL, "%s", buf);

  free(buf);
}

bool bmm_gl2_compile(GLuint *const pshader,
    GLenum const type, char const *const path) {
  GLuint const shader = glCreateShader(type);
  if (shader == 0)
    return false;

  size_t size;
  void *const buf = bmm_io_conts(&size, path);
  if (buf == NULL) {
    glDeleteShader(shader);

    return false;
  }

  GLchar const *str = buf;
  GLint const length = (GLint) size;
  glShaderSource(shader, 1, &str, &length);

  free(buf);

  glCompileShader(shader);

  GLint status;
  glGetShaderiv(shader, GL_COMPILE_STATUS, &status);
  if (status == 0) {
    bmm_gl2_tle(shader, glGetShaderiv, glGetShaderInfoLog);

    glDeleteShader(shader);

    return false;
  }

  *pshader = shader;

  return true;
}

bool bmm_gl2_link(GLuint *const pprogram,
    GLuint const vshader, GLuint const fshader) {
  GLuint const program = glCreateProgram();
  if (program == 0)
    return false;

  glAttachShader(program, vshader);
  glAttachShader(program, fshader);

  glLinkProgram(program);

  GLint status;
  glGetProgramiv(program, GL_LINK_STATUS, &status);
  if (status == 0) {
    bmm_gl2_tle(program, glGetProgramiv, glGetProgramInfoLog);

    glDeleteProgram(program);

    return false;
  }

  *pprogram = program;

  return true;
}

bool bmm_gl2_build(GLuint *const pvshader, GLuint *const pfshader,
    GLuint *const pprogram, char const *const vpath, char const *const fpath) {
  GLuint vshader;
  if (!bmm_gl2_compile(&vshader, GL_VERTEX_SHADER, vpath))
    return false;

  GLuint fshader;
  if (!bmm_gl2_compile(&fshader, GL_FRAGMENT_SHADER, fpath)) {
    glDeleteShader(vshader);

    return false;
  }

  GLuint program;
  if (!bmm_gl2_link(&program, vshader, fshader)) {
    glDeleteShader(fshader);
    glDeleteShader(vshader);

    return false;
  }

  *pvshader = vshader;
  *pfshader = fshader;
  *pprogram = program;

  return true;
}

bool bmm_gl2_buffer(GLuint *const pbuffer,
    GLenum const target, GLenum const usage,
    GLvoid const *const data, GLsizei const size) {
  GLuint buffer;
  glGenBuffers(1, &buffer);

  glBindBuffer(target, buffer);

  glBufferData(target, size, data, usage);
  if (glGetError() != GL_NO_ERROR) {
    glDeleteBuffers(1, &buffer);

    return false;
  }

  *pbuffer = buffer;

  return true;
}
