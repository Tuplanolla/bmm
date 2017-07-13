#ifndef BMM_HIST_H
#define BMM_HIST_H

#include <stddef.h>

#include "ext.h"

// TODO Get rid of useless indirection and enable inlining.
// TODO Restore previous version without weights.
// TODO Give up trying to fix this shit.

/// This structure holds hit counts and cached values.
struct bmm_hist {
  size_t* m;
  double* wm;
  size_t* tmp;
  size_t ndim;
  size_t nlin;
  size_t ncub;
  double a;
  double b;
  size_t nsum;
  double wnsum;
};

/// The call `bmm_hist_free(hist)` releases the memory
/// used by the histogram `hist`.
void bmm_hist_free(struct bmm_hist*);

/// The call `bmm_hist_forget(hist)` clears the histogram `hist`.
__attribute__ ((__nonnull__))
void bmm_hist_forget(struct bmm_hist*);

/// The statement `hist = bmm_hist_alloc(ndim, nlin, a, b)` allocates memory
/// for managing the histogram `hist` of `ndim` dimensions,
/// `nlin` bins per dimension and a range from `a` to `b`.
/// When `hist` is no longer used, `bmm_hist_free` should be called on it.
__attribute__ ((__malloc__))
struct bmm_hist* bmm_hist_alloc(size_t, size_t, double, double);

/// The call `bmm_hist_bin(hist, x)` returns
/// the index of the bin for `x` in the histogram `hist`
/// if `x` is in the valid range.
/// Otherwise `SIZE_MAX` is returned.
/// The time complexity is $O(1)$.
__attribute__ ((__nonnull__, __pure__))
size_t bmm_hist_bin(struct bmm_hist const*, double const*);

/// The call `bmm_hist_unbin(hist, x, i)` writes
/// the center of the bin `i` from the histogram `hist` into `x`.
/// The time complexity is $O(1)$.
__attribute__ ((__nonnull__))
void bmm_hist_unbin(struct bmm_hist const*, double*, size_t);

/// The call `bmm_hist_funbin(hist, x, i)` writes
/// the minimum of the bin `i` from the histogram `hist` into `x`.
/// The time complexity is $O(1)$.
__attribute__ ((__nonnull__))
void bmm_hist_funbin(struct bmm_hist const*, double*, size_t);

/// The call `bmm_hist_cunbin(hist, x, i)` writes
/// the maximum of the bin `i` from the histogram `hist` into `x`.
/// The time complexity is $O(1)$.
__attribute__ ((__nonnull__))
void bmm_hist_cunbin(struct bmm_hist const*, double*, size_t);

/// The call `bmm_hist_accum(hist, x)` adds `x` to the histogram `hist` and
/// returns `true` if `x` is in the valid range.
/// Otherwise `false` is returned.
/// The time complexity is $O(1)$.
__attribute__ ((__nonnull__))
bool bmm_hist_accum(struct bmm_hist*, double const*);

__attribute__ ((__nonnull__))
bool bmm_whist_accum(struct bmm_hist*, double const*, double);

/// The call `bmm_hist_ndim(hist)` returns
/// the number of dimensions in the histogram `hist`.
/// The time complexity is $O(1)$.
__attribute__ ((__nonnull__, __pure__))
size_t bmm_hist_ndim(struct bmm_hist const*);

/// The call `bmm_hist_nsubdiv(hist)` returns
/// the number of bins per dimension in the histogram `hist`.
/// The time complexity is $O(1)$.
__attribute__ ((__nonnull__, __pure__))
size_t bmm_hist_nsubdiv(struct bmm_hist const*);

/// The call `bmm_hist_nbin(hist)` returns
/// the number of bins in the histogram `hist`.
/// The time complexity is $O(1)$.
__attribute__ ((__nonnull__, __pure__))
size_t bmm_hist_nbin(struct bmm_hist const*);

/// The call `bmm_hist_min(hist)` returns
/// the lower limit of the valid range of the histogram `hist`.
/// The time complexity is $O(1)$.
__attribute__ ((__nonnull__, __pure__))
double bmm_hist_min(struct bmm_hist const*);

/// The call `bmm_hist_max(hist)` returns
/// the higher limit of the valid range of the histogram `hist`.
/// The time complexity is $O(1)$.
__attribute__ ((__nonnull__, __pure__))
double bmm_hist_max(struct bmm_hist const*);

/// The call `bmm_hist_length(hist)` returns
/// the length of one dimension of the histogram `hist`.
/// The time complexity is $O(1)$.
__attribute__ ((__nonnull__, __pure__))
double bmm_hist_length(struct bmm_hist const*);

/// The call `bmm_hist_volume(hist)` returns
/// the product of the length of all the dimensions of the histogram `hist`.
/// The time complexity is $O(1)$.
__attribute__ ((__nonnull__, __pure__))
double bmm_hist_volume(struct bmm_hist const*);

/// The call `bmm_hist_hits(hist, i)` returns
/// the hit count for bin `i` in the histogram `hist`.
/// The time complexity is $O(1)$.
__attribute__ ((__nonnull__, __pure__))
size_t bmm_hist_hits(struct bmm_hist const*, size_t);

__attribute__ ((__nonnull__, __pure__))
double bmm_whist_hits(struct bmm_hist const*, size_t);

/// The call `bmm_hist_sumhits(hist)` returns
/// the total hit count in the histogram `hist`.
/// The time complexity is $O(1)$.
__attribute__ ((__nonnull__, __pure__))
size_t bmm_hist_sumhits(struct bmm_hist const*);

__attribute__ ((__nonnull__, __pure__))
double bmm_whist_sumhits(struct bmm_hist const*);

/// The call `bmm_hist_normhits(hist, i)` returns
/// the normalized hit count for bin `i` in the histogram `hist`.
/// The time complexity is $O(1)$.
__attribute__ ((__nonnull__, __pure__))
double bmm_hist_normhits(struct bmm_hist const*, size_t);

#endif
