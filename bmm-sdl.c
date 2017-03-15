#include "opt.h"
#include "ext.h"
#include "sdl.h"
#include "str.h"
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

__attribute__ ((__nonnull__ (1, 2)))
static bool f(char const* const key, char const* const value,
    void* const ptr) {
  struct bmm_sdl_opts* const opts = ptr;

  if (strcmp(key, "width") == 0) {
    if (!bmm_str_strtoz(value, &opts->width))
      return false;
  } else if (strcmp(key, "height") == 0) {
    if (!bmm_str_strtoz(value, &opts->height))
      return false;
  } else if (strcmp(key, "fps") == 0) {
    if (!bmm_str_strtou(value, &opts->fps))
      return false;
  } else if (strcmp(key, "msaa") == 0) {
    if (!bmm_str_strtou(value, &opts->msaa))
      return false;
  }

  return true;
}

__attribute__ ((__nonnull__))
int main(int const argc, char** const argv) {
  struct bmm_sdl_opts opts;
  bmm_sdl_defopts(&opts);
  bmm_opt_parse((char const* const*) &argv[1], (size_t) (argc - 1), f, &opts);

  return bmm_sdl_run(&opts) ? EXIT_SUCCESS : EXIT_FAILURE;
}
