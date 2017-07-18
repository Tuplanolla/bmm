#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "ext.h"
#include "nc.h"
#include "opt.h"
#include "str.h"
#include "tle.h"

__attribute__ ((__nonnull__ (1, 2)))
static bool f(char const *const key, char const *const value,
    void *const ptr) {
  struct bmm_nc_opts *const opts = ptr;

  if (strcmp(key, "conv") == 0) {
    if (strcmp(value, "amber") == 0)
      opts->conv = BMM_NC_CONV_AMBER;
    else
      return false;
  } else if (strcmp(key, "path") == 0) {
    if (strlen(value) < 1)
      return false;

    opts->path = value;
  } else if (strcmp(key, "id") == 0) {
    bool p;
    if (!bmm_str_strtob(&p, value))
      return false;

    opts->i = p;
  } else if (strcmp(key, "radius") == 0) {
    bool p;
    if (!bmm_str_strtob(&p, value))
      return false;

    opts->r = p;
  } else if (strcmp(key, "position") == 0) {
    bool p;
    if (!bmm_str_strtob(&p, value))
      return false;

    opts->q = p;
  } else if (strcmp(key, "velocity") == 0) {
    bool p;
    if (!bmm_str_strtob(&p, value))
      return false;

    opts->v = p;
  } else if (strcmp(key, "force") == 0) {
    bool p;
    if (!bmm_str_strtob(&p, value))
      return false;

    opts->f = p;
  } else
    return false;

  return true;
}

__attribute__ ((__nonnull__))
int main(int const argc, char **const argv) {
  bmm_tle_reset(argv[0]);

  struct bmm_nc_opts opts;
  bmm_nc_opts_def(&opts);

  if (!bmm_opt_parse((char const *const *) &argv[1], (size_t) (argc - 1),
        f, &opts)) {
    bmm_tle_put();

    return EXIT_FAILURE;
  }

  if (!bmm_nc_run_with(&opts)) {
    bmm_tle_put();

    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
