#include <GL/glew.h>

#ifdef FREE
#include <GL/freeglut.h>
#else
#include <GL/glut.h>
#endif

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// TODO Put these into `gl2.h`.

#include "io.h"
#include "tle.h"

__attribute__ ((__nonnull__))
void bmm_gl2_tle(GLuint, PFNGLGETSHADERIVPROC, PFNGLGETSHADERINFOLOGPROC);

__attribute__ ((__nonnull__))
bool bmm_gl2_compile(GLuint*, GLenum, char const*);

__attribute__ ((__nonnull__))
bool bmm_gl2_link(GLuint*, GLuint, GLuint);

__attribute__ ((__nonnull__))
bool bmm_gl2_build(GLuint*, GLuint*, GLuint*, char const*, char const*);

void bmm_gl2_tle(GLuint const object,
    PFNGLGETSHADERIVPROC const iv,
    PFNGLGETSHADERINFOLOGPROC const infoLog) {
  GLint length;
  iv(object, GL_INFO_LOG_LENGTH, &length);
  if (length <= 0)
    return;

  size_t const size = (size_t) length;

  char* const buf = malloc(size);
  if (buf == NULL) {
    BMM_TLE_STDS();

    return;
  }

  infoLog(object, (GLsizei) size, NULL, (GLchar*) buf);
  BMM_TLE_EXTS(BMM_TLE_NUM_GL, "%s", buf);

  free(buf);
}

bool bmm_gl2_compile(GLuint* const pshader,
    GLenum const type, char const* const path) {
  GLuint const shader = glCreateShader(type);
  if (shader == 0)
    return false;

  size_t size;
  void* const buf = bmm_io_conts(&size, path);
  if (buf == NULL) {
    glDeleteShader(shader);

    return false;
  }

  GLchar const* const str = buf;
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

bool bmm_gl2_link(GLuint* const pprogram,
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

bool bmm_gl2_build(GLuint* const pprogram,
    GLuint* const pvshader, GLuint* const pfshader,
    char const* const vpath, char const* const fpath) {
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

  *pprogram = program;
  *pvshader = vshader;
  *pfshader = fshader;

  return true;
}

// Thank you, `www.opengl-tutorial.org` and `learnopengl.com`.

// Dirty.
static struct {
  GLuint vertex_buffer, element_buffer;

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
static bool make_resources(void) {
  g_resources.vertex_buffer = make_buffer(
      GL_ARRAY_BUFFER,
      g_vertex_buffer_data,
      sizeof g_vertex_buffer_data);
  g_resources.element_buffer = make_buffer(
      GL_ELEMENT_ARRAY_BUFFER,
      g_element_buffer_data,
      sizeof g_element_buffer_data);

  if (!bmm_gl2_build(&g_resources.program,
      &g_resources.vertex_shader,
      &g_resources.fragment_shader,
      "glut.vs.glsl", "glut.fs.glsl"))
    return false;

  g_resources.uniforms.fade_factor
    = glGetUniformLocation(g_resources.program, "fade_factor");

  g_resources.attributes.position
    = glGetAttribLocation(g_resources.program, "position");

  return true;
}

// Dirty.
static void render(void) {

#if 1

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

#else

  // Clear the screen
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  double delta = 0.1;

  computeMatricesFromInputs();
  glm::mat4 ProjectionMatrix = getProjectionMatrix();
  glm::mat4 ViewMatrix = getViewMatrix();

  // We will need the camera's position in order to sort the particles
  // w.r.t the camera's distance.
  // There should be a getCameraPosition() function in common/controls.cpp, 
  // but this works too.
  glm::vec3 CameraPosition(glm::inverse(ViewMatrix)[3]);

  glm::mat4 ViewProjectionMatrix = ProjectionMatrix * ViewMatrix;

  // Generate 10 new particule each millisecond,
  // but limit this to 16 ms (60 fps), or if you have 1 long frame (1sec),
  // newparticles will be huge and the next frame even longer.
  int newparticles = (int)(delta*10000.0);
  if (newparticles > (int)(0.016f*10000.0))
    newparticles = (int)(0.016f*10000.0);

  for(int i=0; i<newparticles; i++){
    int particleIndex = FindUnusedParticle();
    ParticlesContainer[particleIndex].life = 5.0f; // This particle will live 5 seconds.
    ParticlesContainer[particleIndex].pos = glm::vec3(0,0,-20.0f);

    float spread = 1.5f;
    glm::vec3 maindir = glm::vec3(0.0f, 10.0f, 0.0f);
    // Very bad way to generate a random direction; 
    // See for instance http://stackoverflow.com/questions/5408276/python-uniform-spherical-distribution instead,
    // combined with some user-controlled parameters (main direction, spread, etc)
    glm::vec3 randomdir = glm::vec3(
        (rand()%2000 - 1000.0f)/1000.0f,
        (rand()%2000 - 1000.0f)/1000.0f,
        (rand()%2000 - 1000.0f)/1000.0f
        );

    ParticlesContainer[particleIndex].speed = maindir + randomdir*spread;

    // Very bad way to generate a random color
    ParticlesContainer[particleIndex].r = rand() % 256;
    ParticlesContainer[particleIndex].g = rand() % 256;
    ParticlesContainer[particleIndex].b = rand() % 256;
    ParticlesContainer[particleIndex].a = (rand() % 256) / 3;

    ParticlesContainer[particleIndex].size = (rand()%1000)/2000.0f + 0.1f;

  }

  // Simulate all particles
  int ParticlesCount = 0;
  for(int i=0; i<MaxParticles; i++){

    Particle& p = ParticlesContainer[i]; // shortcut

    if(p.life > 0.0f){

      // Decrease life
      p.life -= delta;
      if (p.life > 0.0f){

        // Simulate simple physics : gravity only, no collisions
        p.speed += glm::vec3(0.0f,-9.81f, 0.0f) * (float)delta * 0.5f;
        p.pos += p.speed * (float)delta;
        p.cameradistance = glm::length2( p.pos - CameraPosition );
        //ParticlesContainer[i].pos += glm::vec3(0.0f,10.0f, 0.0f) * (float)delta;

        // Fill the GPU buffer
        g_particule_position_size_data[4*ParticlesCount+0] = p.pos.x;
        g_particule_position_size_data[4*ParticlesCount+1] = p.pos.y;
        g_particule_position_size_data[4*ParticlesCount+2] = p.pos.z;

        g_particule_position_size_data[4*ParticlesCount+3] = p.size;

        g_particule_color_data[4*ParticlesCount+0] = p.r;
        g_particule_color_data[4*ParticlesCount+1] = p.g;
        g_particule_color_data[4*ParticlesCount+2] = p.b;
        g_particule_color_data[4*ParticlesCount+3] = p.a;

      }else{
        // Particles that just died will be put at the end of the buffer in SortParticles();
        p.cameradistance = -1.0f;
      }

      ParticlesCount++;

    }
  }

  SortParticles();
  //printf("%d ",ParticlesCount);

  // Update the buffers that OpenGL uses for rendering.
  // There are much more sophisticated means to stream data from the CPU to the GPU, 
  // but this is outside the scope of this tutorial.
  // http://www.opengl.org/wiki/Buffer_Object_Streaming

  glBindBuffer(GL_ARRAY_BUFFER, particles_position_buffer);
  glBufferData(GL_ARRAY_BUFFER, MaxParticles * 4 * sizeof(GLfloat), NULL, GL_STREAM_DRAW); // Buffer orphaning, a common way to improve streaming perf. See above link for details.
  glBufferSubData(GL_ARRAY_BUFFER, 0, ParticlesCount * sizeof(GLfloat) * 4, g_particule_position_size_data);

  glBindBuffer(GL_ARRAY_BUFFER, particles_color_buffer);
  glBufferData(GL_ARRAY_BUFFER, MaxParticles * 4 * sizeof(GLubyte), NULL, GL_STREAM_DRAW); // Buffer orphaning, a common way to improve streaming perf. See above link for details.
  glBufferSubData(GL_ARRAY_BUFFER, 0, ParticlesCount * sizeof(GLubyte) * 4, g_particule_color_data);


  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  // Use our shader
  glUseProgram(programID);

  // Same as the billboards tutorial
  glUniform3f(CameraRight_worldspace_ID, ViewMatrix[0][0], ViewMatrix[1][0], ViewMatrix[2][0]);
  glUniform3f(CameraUp_worldspace_ID   , ViewMatrix[0][1], ViewMatrix[1][1], ViewMatrix[2][1]);

  glUniformMatrix4fv(ViewProjMatrixID, 1, GL_FALSE, &ViewProjectionMatrix[0][0]);

  // 1rst attribute buffer : vertices
  glEnableVertexAttribArray(0);
  glBindBuffer(GL_ARRAY_BUFFER, billboard_vertex_buffer);
  glVertexAttribPointer(
      0,                  // attribute. No particular reason for 0, but must match the layout in the shader.
      3,                  // size
      GL_FLOAT,           // type
      GL_FALSE,           // normalized?
      0,                  // stride
      (void*)0            // array buffer offset
      );

  // 2nd attribute buffer : positions of particles' centers
  glEnableVertexAttribArray(1);
  glBindBuffer(GL_ARRAY_BUFFER, particles_position_buffer);
  glVertexAttribPointer(
      1,                                // attribute. No particular reason for 1, but must match the layout in the shader.
      4,                                // size : x + y + z + size => 4
      GL_FLOAT,                         // type
      GL_FALSE,                         // normalized?
      0,                                // stride
      (void*)0                          // array buffer offset
      );

  // 3rd attribute buffer : particles' colors
  glEnableVertexAttribArray(2);
  glBindBuffer(GL_ARRAY_BUFFER, particles_color_buffer);
  glVertexAttribPointer(
      2,                                // attribute. No particular reason for 1, but must match the layout in the shader.
      4,                                // size : r + g + b + a => 4
      GL_UNSIGNED_BYTE,                 // type
      GL_TRUE,                          // normalized?    *** YES, this means that the unsigned char[4] will be accessible with a vec4 (floats) in the shader ***
      0,                                // stride
      (void*)0                          // array buffer offset
      );

  // These functions are specific to glDrawArrays*Instanced*.
  // The first parameter is the attribute buffer we're talking about.
  // The second parameter is the "rate at which generic vertex attributes advance when rendering multiple instances"
  // http://www.opengl.org/sdk/docs/man/xhtml/glVertexAttribDivisor.xml
  glVertexAttribDivisor(0, 0); // particles vertices : always reuse the same 4 vertices -> 0
  glVertexAttribDivisor(1, 1); // positions : one per quad (its center)                 -> 1
  glVertexAttribDivisor(2, 1); // color : one per quad                                  -> 1

  // Draw the particules !
  // This draws many times a small triangle_strip (which looks like a quad).
  // This is equivalent to :
  // for(i in ParticlesCount) : glDrawArrays(GL_TRIANGLE_STRIP, 0, 4), 
  // but faster.
  glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, ParticlesCount);

  glDisableVertexAttribArray(0);
  glDisableVertexAttribArray(1);
  glDisableVertexAttribArray(2);

#endif

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

#if 0

  GLuint VertexArrayID;
  glGenVertexArrays(1, &VertexArrayID);
  glBindVertexArray(VertexArrayID);

  // Create and compile our GLSL program from the shaders
  GLuint programID = LoadShaders( "Particle.vertexshader", "Particle.fragmentshader" );

  // Vertex shader
  GLuint CameraRight_worldspace_ID  = glGetUniformLocation(programID, "CameraRight_worldspace");
  GLuint CameraUp_worldspace_ID  = glGetUniformLocation(programID, "CameraUp_worldspace");
  GLuint ViewProjMatrixID = glGetUniformLocation(programID, "VP");

  static GLfloat* g_particule_position_size_data = new GLfloat[MaxParticles * 4];
  static GLubyte* g_particule_color_data         = new GLubyte[MaxParticles * 4];

  for(int i=0; i<MaxParticles; i++){
    ParticlesContainer[i].life = -1.0f;
    ParticlesContainer[i].cameradistance = -1.0f;
  }

  // The VBO containing the 4 vertices of the particles.
  // Thanks to instancing, they will be shared by all particles.
  static const GLfloat g_vertex_buffer_data[] = { 
    -0.5f, -0.5f, 0.0f,
    0.5f, -0.5f, 0.0f,
    -0.5f,  0.5f, 0.0f,
    0.5f,  0.5f, 0.0f,
  };
  GLuint billboard_vertex_buffer;
  glGenBuffers(1, &billboard_vertex_buffer);
  glBindBuffer(GL_ARRAY_BUFFER, billboard_vertex_buffer);
  glBufferData(GL_ARRAY_BUFFER, sizeof(g_vertex_buffer_data), g_vertex_buffer_data, GL_STATIC_DRAW);

  // The VBO containing the positions and sizes of the particles
  GLuint particles_position_buffer;
  glGenBuffers(1, &particles_position_buffer);
  glBindBuffer(GL_ARRAY_BUFFER, particles_position_buffer);
  // Initialize with empty (NULL) buffer : it will be updated later, each frame.
  glBufferData(GL_ARRAY_BUFFER, MaxParticles * 4 * sizeof(GLfloat), NULL, GL_STREAM_DRAW);

  // The VBO containing the colors of the particles
  GLuint particles_color_buffer;
  glGenBuffers(1, &particles_color_buffer);
  glBindBuffer(GL_ARRAY_BUFFER, particles_color_buffer);
  // Initialize with empty (NULL) buffer : it will be updated later, each frame.
  glBufferData(GL_ARRAY_BUFFER, MaxParticles * 4 * sizeof(GLubyte), NULL, GL_STREAM_DRAW);

  delete[] g_particule_position_size_data;

  // Cleanup VBO and shader
  glDeleteBuffers(1, &particles_color_buffer);
  glDeleteBuffers(1, &particles_position_buffer);
  glDeleteBuffers(1, &billboard_vertex_buffer);
  glDeleteProgram(programID);
  glDeleteVertexArrays(1, &VertexArrayID);

#endif

  glutMainLoop();

  // TODO Nope.
  fputs(bmm_tle_msg(), stderr);

#ifdef FREE
  before_exit();

  return EXIT_SUCCESS;
#else
  return EXIT_FAILURE;
#endif
}
