#include "glut.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "ext.h"
#include "opt.h"
#include "str.h"
#include "tle.h"

__attribute__ ((__nonnull__ (1, 2)))
static bool f(char const* const key, char const* const value,
    void* const ptr) {
  struct bmm_glut_opts* const opts = ptr;

  if (strcmp(key, "vpath") == 0) {
    if (strlen(value) < 1)
      return false;

    opts->vpath = value;
  } else if (strcmp(key, "fpath") == 0) {
    if (strlen(value) < 1)
      return false;

    opts->fpath = value;
  } else
    return false;

  return true;
}

__attribute__ ((__nonnull__))
int main(int const argc, char** const argv) {
  bmm_tle_reset(argv[0]);

  struct bmm_glut_opts opts;
  bmm_glut_opts_def(&opts);

  if (!bmm_opt_parse((char const* const*) &argv[1], (size_t) (argc - 1),
        f, &opts)) {
    bmm_tle_put();

    return EXIT_FAILURE;
  }

  if (!bmm_glut_run_with(&opts)) {
    bmm_tle_put();

    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
