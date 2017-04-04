#include "err.h"
#include "ext.h"
#include "filter.h"
#include "opt.h"
#include <stdlib.h>
#include <string.h>

// TODO Rethink this user interface.

__attribute__ ((__nonnull__ (1, 2)))
static bool f(char const* const key, char const* const value,
    void* const ptr) {
  struct bmm_filter_opts* const opts = ptr;

  if (strcmp(key, "basis") == 0) {
    size_t i;
    if (!bmm_str_strtoz(&i, value))
      return false;
    for (size_t imsg = 0; imsg < BMM_MSG_MAX; ++imsg)
      opts->mask[imsg] = i != 0;
  } else if (strcmp(key, "with") == 0) {
    size_t i;
    if (!bmm_str_strtoz(&i, value))
      return false;
    opts->mask[i] = true;
  } else if (strcmp(key, "without") == 0) {
    size_t i;
    if (!bmm_str_strtoz(&i, value))
      return false;
    opts->mask[i] = false;
  } else
    return false;

  return true;
}

__attribute__ ((__nonnull__))
int main(int const argc, char** const argv) {
  bmm_err_reset();

  struct bmm_filter_opts opts;
  bmm_filter_defopts(&opts);
  bmm_opt_parse((char const* const*) &argv[1], (size_t) (argc - 1), f, &opts);

  return bmm_filter_run(&opts) ? EXIT_SUCCESS : EXIT_FAILURE;
}
