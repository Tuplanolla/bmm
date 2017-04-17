#version 330 core

uniform float fade_factor;

void main(void) {
  gl_FragColor = vec4(fade_factor);
}
