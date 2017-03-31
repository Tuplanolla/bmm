#include "dem.h"
#include "err.h"
#include "filter.h"
#include "msg.h"
#include <stdio.h>

void bmm_filter_defopts(struct bmm_filter_opts* const opts) {
  opts->dummy = 0;
}

void bmm_filter_def(struct bmm_filter* const filter,
    struct bmm_filter_opts const* const opts) {
  struct bmm_dem_opts defopts;
  bmm_dem_defopts(&defopts);
  bmm_dem_def(&filter->dem, &defopts);

  filter->opts = *opts;
}

static bool f(struct bmm_msg_head const* const head,
    void* const ptr) {
  bool* const eof = ptr;
  *eof = false;

  // TODO Think about this.
  // return head->type == BMM_MSG_EKINE;
  return true;
}

bool bmm_filter_run_with(struct bmm_filter* const filter) {
  for ever {
    bool eof = true;

    struct bmm_msg_head head;
    if (!bmm_msg_get(&head, &filter->dem, f, &eof))
      return false;

    if (eof)
      return true;

    bmm_msg_put(&head, &filter->dem);

    fprintf(stderr, "%g %g %g %g\n",
        filter->dem.opts.tstep * (double) filter->dem.istep,
        filter->dem.est.ekinetic,
        filter->dem.est.pvector, filter->dem.est.pscalar);
  }

  return true;
}

bool bmm_filter_run(struct bmm_filter_opts const* const opts) {
  struct bmm_filter* const filter = malloc(sizeof *filter);
  if (filter == NULL) {
    BMM_ERR_WARN(malloc);

    return false;
  }

  bmm_filter_def(filter, opts);
  bool const result = bmm_filter_run_with(filter);

  free(filter);

  return result;
}
