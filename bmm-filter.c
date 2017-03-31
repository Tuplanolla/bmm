#include "err.h"
#include "ext.h"
#include "filter.h"
#include <stdlib.h>

__attribute__ ((__nonnull__))
int main(int const argc, char** const argv) {
  bmm_err_reset();

  struct bmm_filter_opts opts;
  bmm_filter_defopts(&opts);

  return bmm_filter_run(&opts) ? EXIT_SUCCESS : EXIT_FAILURE;
}
