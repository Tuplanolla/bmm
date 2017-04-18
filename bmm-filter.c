#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "conf.h"
#include "ext.h"
#include "filter.h"
#include "opt.h"
#include "str.h"
#include "tle.h"

__attribute__ ((__nonnull__ (1, 2)))
static bool f(char const* const key, char const* const value,
    void* const ptr) {
  struct bmm_filter_opts* const opts = ptr;

  if (strcmp(key, "mode") == 0) {
    if (strcmp(value, "blacklist") == 0)
      for (size_t imsg = 0; imsg < BMM_MSG_MAX; ++imsg)
        opts->mask[imsg] = true;
    else if (strcmp(value, "whitelist") == 0)
      for (size_t imsg = 0; imsg < BMM_MSG_MAX; ++imsg)
        opts->mask[imsg] = false;
    else
      return false;
  } else if (strcmp(key, "pass") == 0) {
    size_t i;
    if (!(bmm_str_strtoz(&i, value) && i < BMM_MSG_MAX))
      return false;
    opts->mask[i] = true;
  } else if (strcmp(key, "stop") == 0) {
    size_t i;
    if (!(bmm_str_strtoz(&i, value) && i < BMM_MSG_MAX))
      return false;
    opts->mask[i] = false;
  } else if (strcmp(key, "verbose") == 0) {
    bool p;
    if (!bmm_str_strtob(&p, value))
      return false;
    opts->verbose = p;
  } else
    return false;

  return true;
}

__attribute__ ((__nonnull__))
int main(int const argc, char** const argv) {
  bmm_tle_reset(argv[0]);

  struct bmm_filter_opts opts;
  bmm_filter_opts_def(&opts);

  if (!bmm_opt_parse((char const* const*) &argv[1], (size_t) (argc - 1),
        f, &opts)) {
    bmm_tle_put();

    return EXIT_FAILURE;
  }

  if (!bmm_filter_run_with(&opts)) {
    bmm_tle_put();

    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
