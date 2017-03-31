// Message filter.
#ifndef BMM_FILTER_H
#define BMM_FILTER_H

#include "dem.h"
#include "ext.h"
#include <limits.h>
#include <stdbool.h>

struct bmm_filter_opts {
  int dummy;
};

struct bmm_filter {
  struct bmm_filter_opts opts;
  struct bmm_dem dem;
};

__attribute__ ((__nonnull__))
void bmm_filter_defopts(struct bmm_filter_opts*);

__attribute__ ((__nonnull__))
void bmm_filter_def(struct bmm_filter*, struct bmm_filter_opts const*);

__attribute__ ((__nonnull__))
bool bmm_filter_run_with(struct bmm_filter*);

__attribute__ ((__nonnull__))
bool bmm_filter_run(struct bmm_filter_opts const*);

#endif
