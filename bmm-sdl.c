#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "err.h"
#include "ext.h"
#include "opt.h"
#include "sdl.h"
#include "str.h"

__attribute__ ((__nonnull__ (1, 2)))
static bool f(char const* const key, char const* const value,
    void* const ptr) {
  struct bmm_sdl_opts* const opts = ptr;

  if (strcmp(key, "width") == 0) {
    if (!bmm_str_strtou(&opts->width, value))
      return false;
  } else if (strcmp(key, "height") == 0) {
    if (!bmm_str_strtou(&opts->height, value))
      return false;
  } else if (strcmp(key, "fps") == 0) {
    if (!bmm_str_strtou(&opts->fps, value))
      return false;
  } else if (strcmp(key, "ms") == 0) {
    if (!bmm_str_strtou(&opts->ms, value))
      return false;
  } else if (strcmp(key, "zoomfac") == 0) {
    if (!bmm_str_strtod(&opts->zoomfac, value))
      return false;
  } else
    return false;

  return true;
}

__attribute__ ((__nonnull__))
int main(int const argc, char** const argv) {
  bmm_tle_reset(argv[0]);

  struct bmm_sdl_opts opts;
  bmm_sdl_opts_def(&opts);

  if (!bmm_opt_parse((char const* const*) &argv[1], (size_t) (argc - 1),
        f, &opts))
    return EXIT_FAILURE;

  if (!bmm_sdl_run_with(&opts))
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
