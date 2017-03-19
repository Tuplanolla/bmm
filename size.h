// Index space operations.
#ifndef BMM_SIZE_H
#define BMM_SIZE_H

#include "ext.h"
#include <limits.h>
#include <stddef.h>

// This structure holds the quotient and remainder of a division
// in unspecified order.
typedef struct {
  size_t quot;
  size_t rem;
} bmm_size_div_t;

// The call `m = bmm_size_div(n, k)` solves
// the division equation `m.quot * k + m.rem == n` for `m`,
// where `m.quot` is the quotient and `m.rem` is the remainder
// of the expression `n / k`.
// This is analogous to `div`.
__attribute__ ((__const__, __pure__))
inline bmm_size_div_t bmm_size_div(size_t const n, size_t const k) {
  bmm_size_div_t const qr = {
    .quot = n / k,
    .rem = n % k
  };

  return qr;
}

// The call `bmm_size_cmp(n, k)` returns
//
// * `-1` if `n < k`,
// * `1` if `n > k` and
// * `0` otherwise.
//
// This is analogous to `bmm_fp_cmp`.
__attribute__ ((__const__, __pure__))
inline int bmm_size_cmp(size_t const n, size_t const k) {
  return n < k ? -1 : n > k ? 1 : 0;
}

// The call `bmm_size_min(n, k)` returns the lesser of `n` and `k`.
// This is analogous to `fmin`.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_min(size_t const n, size_t const k) {
  return n < k ? n : k;
}

// The call `bmm_size_max(n, k)` returns the greater of `n` and `k`.
// This is analogous to `fmax`.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_max(size_t const n, size_t const k) {
  return n > k ? n : k;
}

// The call `bmm_size_pow(n, k)` returns `n` raised to the power of `k`.
// This is analogous to `pow`.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_pow(size_t const n, size_t const k) {
  size_t m = 1;

  for (size_t i = 0; i < k; ++i)
    m *= n;

  return m;
}

// The call `bmm_size_identity(n)` returns `n`.
// This is analogous to `bmm_fp_identity`.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_identity(size_t const n) {
  return n;
}

// The call `bmm_size_constant(n, k)` returns `n`.
// This is analogous to `bmm_fp_constant`.
__attribute__ ((__const__, __pure__))
inline size_t bmm_fp_constant(size_t const x,
    __attribute__ ((__unused__)) size_t const y) {
  return x;
}

// The call `bmm_size_zero(n)` returns `0`.
// This is analogous to `bmm_fp_zero`.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_zero(__attribute__ ((__unused__)) size_t const n) {
  return 0;
}

// The call `bmm_size_one(n)` returns `1`.
// This is analogous to `bmm_fp_one`.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_one(__attribute__ ((__unused__)) size_t const n) {
  return 1;
}

// The call `bmm_size_midpoint(n, k)`
// returns the arithmetic mean of `n` and `k`.
// This is analogous to `bmm_fp_midpoint`.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_midpoint(size_t const n, size_t const k) {
  // Note that `(n + k) / 2` could wrap and
  // `n < k ? n + (k - n) / 2 : k + (n - k) / 2` could have bad performance.
  return n / 2 + k / 2 + (n % 2 + k % 2) / 2;
}

// The call `bmm_size_sq(n)` returns `n` squared.
// This is analogous to `bmm_fp_sq`.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_sq(size_t const n) {
  return n * n;
}

// The call `bmm_size_firt(n, k)` returns the floor of the `k`th root of `n`.
// This is analogous to `bmm_fp_rt`.
// Note that the result may be wrong for large arguments.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_firt(size_t const n, size_t const k) {
  if (n <= 1)
    return n;
  else {
    size_t const p = k - 1;
    size_t r = n + 1;
    size_t m = n;

    while (m < r) {
      r = m;
      m = (p * r + n / bmm_size_pow(r, p)) / k;
    }

    return r;
  }
}

// The call `bmm_size_cirt(n, k)` returns the ceiling of the `k`th root of `n`.
// This is analogous to `bmm_fp_rt`.
// Note that the result may be wrong for large arguments.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_cirt(size_t const n, size_t const k) {
  return n <= 1 ? n : bmm_size_firt(n - 1, k) + 1;
}

// The call `bmm_size_filog(n, k)` returns the floor
// of the base `k` logarithm of `n`.
// This is analogous to `bmm_fp_log`.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_filog(size_t n, size_t const k) {
  dynamic_assert(n <= 0, "invalid argument");
  dynamic_assert(k <= 1, "invalid base");

  size_t m = 0;

  while (n >= k) {
    n /= k;
    ++m;
  }

  return m;
}

// The call `bmm_size_cilog(n, k)` returns the ceiling
// of the base `k` logarithm of `n`.
// This is analogous to `bmm_fp_log`.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_cilog(size_t const n, size_t const k) {
  dynamic_assert(n <= 0, "invalid argument");
  dynamic_assert(k <= 1, "invalid base");

  return n <= 1 ? 0 : bmm_size_filog(n - 1, k) + 1;
}

// The call `bmm_size_uclamp(n, k)` returns
//
// * `n` if `0 <= n < k` and
// * `k - 1` if `n >= k`.
//
// This is analogous to `bmm_fp_uclamp`.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_uclamp(size_t const n, size_t const k) {
  return n >= k ? k - 1 : n;
}

// The call `m = bmm_size_uwrap(n, k)` solves
// the periodic equation `m == n - p * k` for `m`,
// where `0 <= m < k` and `p` is some integer.
// This is analogous to `bmm_fp_uwrap`.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_uwrap(size_t const n, size_t const k) {
  return n % k;
}

// The call `bmm_size_uwrap_inc(n, k)` is equivalent
// to `bmm_size_uwrap(n + 1, k)` without wrapping.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_uwrap_inc(size_t const n, size_t const k) {
  return n == k - 1 ? 0 : n + 1;
}

// The call `bmm_size_uwrap_dec(n, k)` is equivalent
// to `bmm_size_uwrap(n - 1, k)` without wrapping.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_uwrap_dec(size_t const n, size_t const k) {
  return n == 0 ? k - 1 : n - 1;
}

#endif
