// Message filter.
#ifndef BMM_FILTER_H
#define BMM_FILTER_H

#include "conf.h"
#include "dem.h"
#include "ext.h"
#include "msg.h"
#include <stdbool.h>

struct bmm_filter_opts {
  bool mask[BMM_MSG_MAX];
};

struct bmm_filter {
  struct bmm_filter_opts opts;
  size_t passed;
  size_t stopped;
};

// TODO Mention whitelist by default.
__attribute__ ((__nonnull__))
void bmm_filter_opts_def(struct bmm_filter_opts*);

__attribute__ ((__nonnull__))
void bmm_filter_def(struct bmm_filter*, struct bmm_filter_opts const*);

__attribute__ ((__nonnull__))
bool bmm_filter_run(struct bmm_filter*);

__attribute__ ((__nonnull__))
bool bmm_filter_run_with(struct bmm_filter_opts const*);

#endif
