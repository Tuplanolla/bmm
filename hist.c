#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "common.h"
#include "ext.h"
#include "fp.h"
#include "hist.h"

void bmm_hist_free(struct bmm_hist *const hist) {
  if (hist != NULL) {
    free(hist->tmp);

    free(hist->m);
    free(hist->wm);
  }

  free(hist);
}

void bmm_hist_forget(struct bmm_hist *const hist) {
  for (size_t i = 0; i < hist->ncub; ++i)
    hist->m[i] = 0;

  for (size_t i = 0; i < hist->ncub; ++i)
    hist->wm[i] = 0.0;

  hist->nsum = 0;

  hist->wnsum = 0.0;
}

struct bmm_hist *bmm_hist_alloc(size_t const ndim, size_t const nlin,
    double const a, double const b) {
  bool p = true;

  struct bmm_hist *const hist = malloc(sizeof *hist);
  if (hist == NULL)
    p = false;
  else {
    hist->ndim = ndim;
    hist->nlin = nlin;
    hist->ncub = bmm_size_pow(nlin, ndim);
    hist->a = a;
    hist->b = b;

    hist->m = malloc(hist->ncub * sizeof *hist->m);
    if (hist->m == NULL)
      p = false;

    hist->wm = malloc(hist->ncub * sizeof *hist->wm);
    if (hist->wm == NULL)
      p = false;

    hist->tmp = malloc(ndim * sizeof *hist->tmp);
    if (hist->tmp == NULL)
      p = false;

    if (p)
      bmm_hist_forget(hist);
  }

  if (p)
    return hist;
  else {
    bmm_hist_free(hist);

    return NULL;
  }
}

size_t bmm_hist_bin(struct bmm_hist const *const hist, double const *const x) {
  for (size_t i = 0; i < hist->ndim; ++i)
    if (x[i] >= hist->a && x[i] < hist->b)
      // TODO Very bad!
      hist->tmp[i] = bmm_fp_iclerp(x[i], hist->a, hist->b,
          1, hist->nlin);
    else
      return SIZE_MAX;

  return inst(bmm_unhc, size_t)(hist->tmp, hist->ndim, hist->nlin);
}

__attribute__ ((__nonnull__))
static void unbin(struct bmm_hist const *const hist, double *const x,
    size_t const j, double const h) {
  inst(bmm_hc, size_t)(hist->tmp, j, hist->ndim, hist->nlin);

  for (size_t i = 0; i < hist->ndim; ++i)
    x[i] = bmm_fp_lerp((double) hist->tmp[i] + h,
        0.0, (double) hist->nlin, hist->a, hist->b);
}

void bmm_hist_unbin(struct bmm_hist const *const hist, double *const x,
    size_t const j) {
  unbin(hist, x, j, 0.5);
}

void bmm_hist_funbin(struct bmm_hist const *const hist, double *const x,
    size_t const j) {
  unbin(hist, x, j, 0.0);
}

void bmm_hist_cunbin(struct bmm_hist const *const hist, double *const x,
    size_t const j) {
  unbin(hist, x, j, 1.0);
}

bool bmm_hist_accum(struct bmm_hist *const hist, double const *const x) {
  size_t const i = bmm_hist_bin(hist, x);

  if (i == SIZE_MAX)
    return false;
  else {
    ++hist->m[i];

    ++hist->nsum;

    return true;
  }
}

bool bmm_whist_accum(struct bmm_hist *const hist, double const *const x,
    double const w) {
  size_t const i = bmm_hist_bin(hist, x);

  if (i == SIZE_MAX)
    return false;
  else {
    hist->wm[i] += w;

    hist->wnsum += w;

    return true;
  }
}

size_t bmm_hist_ndim(struct bmm_hist const *const hist) {
  return hist->ndim;
}

size_t bmm_hist_nsubdiv(struct bmm_hist const *const hist) {
  return hist->nlin;
}

size_t bmm_hist_nbin(struct bmm_hist const *const hist) {
  return hist->ncub;
}

double bmm_hist_min(struct bmm_hist const *const hist) {
  return hist->a;
}

double bmm_hist_max(struct bmm_hist const *const hist) {
  return hist->b;
}

double bmm_hist_length(struct bmm_hist const *const hist) {
  return hist->b - hist->a;
}

double bmm_hist_volume(struct bmm_hist const *const hist) {
  return bmm_fp_pow(hist->b - hist->a, hist->ndim);
}

size_t bmm_hist_hits(struct bmm_hist const *const hist, size_t const i) {
  return hist->m[i];
}

double bmm_whist_hits(struct bmm_hist const *const hist, size_t const i) {
  return hist->wm[i];
}

size_t bmm_hist_sumhits(struct bmm_hist const *const hist) {
  return hist->nsum;
}

double bmm_whist_sumhits(struct bmm_hist const *const hist) {
  return hist->wnsum;
}

double bmm_hist_normhits(struct bmm_hist const *const hist, size_t const i) {
  return hist->nsum == 0 ? 0.0 : (double) hist->m[i] /
    (bmm_hist_volume(hist) * ((double) hist->nsum / (double) hist->ncub));
}
