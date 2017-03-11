#include "clopts.h"
#include "exts.h"
#include "strs.h"
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

struct bmm_opts {
  size_t npart;
  size_t nstep;
};

static bool f(char const* const key, char const* const value, void* const p) {
  struct bmm_opts* const opts = p;

  if (strcmp(key, "npart") == 0) {
    if (!bmm_strtoz(value, &opts->npart))
      return false;
  }

  return true;
}

__attribute__ ((__nonnull__))
int main(int const argc, char** const argv) {
  struct bmm_opts opts;
  opts.npart = 0;
  opts.nstep = 0;

  bmm_clopts((char const* const*) &argv[1], (size_t) (argc - 1), f, &opts);

  return EXIT_SUCCESS;
}
