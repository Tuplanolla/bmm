#ifndef BMM_FILTER_H
/// Message filter.
#define BMM_FILTER_H

#include <stdbool.h>
#include <stddef.h>

#include "conf.h"
#include "ext.h"
#include "io.h"

/// This structure contains filter options such as the whitelist.
struct bmm_filter_opts {
  bool verbose;
  bool mask[BMM_MMSG];
};

/// This structure holds some filter statistics.
struct bmm_filter {
  struct bmm_filter_opts opts;
  size_t passed;
  size_t stopped;
};

/// The call `bmm_filter_opts_def(opts)`
/// writes the default filter options into `opts`.
/// All messages are stopped by default.
__attribute__ ((__nonnull__))
void bmm_filter_opts_def(struct bmm_filter_opts*);

/// The call `bmm_filter_def(filter, opts)`
/// writes the default filter state into `filter`
/// with the filter options `opts`.
__attribute__ ((__nonnull__))
void bmm_filter_def(struct bmm_filter*, struct bmm_filter_opts const*);

/// The call `bmm_filter_step(filter)`
/// processes one incoming message with the filter state `filter`.
__attribute__ ((__nonnull__))
enum bmm_io_read bmm_filter_step(struct bmm_filter*);

/// The call `bmm_filter_report(filter)`
/// prints informal diagnostics for the filter state `filter`.
__attribute__ ((__nonnull__))
bool bmm_filter_report(struct bmm_filter const*);

/// The call `bmm_filter_run(filter)`
/// processes all incoming messages and
/// handles signals with the filter state `filter`.
__attribute__ ((__nonnull__))
bool bmm_filter_run(struct bmm_filter*);

/// The call `bmm_filter_run_with(opts)`
/// processes all incoming messages and
/// handles signals with the filter options `opts`.
__attribute__ ((__nonnull__))
bool bmm_filter_run_with(struct bmm_filter_opts const*);

#endif
