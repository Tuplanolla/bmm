#ifndef BMM_SIZE_H
/// Index space operations.
#define BMM_SIZE_H

#include <stddef.h>
#include <stdint.h>

#include "ext.h"

/// This structure holds the quotient and remainder of a division
/// in unspecified order.
typedef struct {
  size_t quot;
  size_t rem;
} bmm_size_div_t;

/// The call `m = bmm_size_div(n, k)`
/// solves the division equation `m.quot * k + m.rem == n` for `m`,
/// where `m.quot` is the quotient and `m.rem` is the remainder
/// of the expression `n / k`.
/// This is analogous to `div` or `bmm_fp_div`.
__attribute__ ((__const__, __pure__))
inline bmm_size_div_t bmm_size_div(size_t const n, size_t const k) {
  bmm_size_div_t const qr = {
    .quot = n / k,
    .rem = n % k
  };

  return qr;
}

/// The call `bmm_size_cmp(n, k)` returns
///
/// * `-1` if `n < k`,
/// * `1` if `n > k` and
/// * `0` otherwise.
///
/// This is analogous to `bmm_fp_cmp`.
__attribute__ ((__const__, __pure__))
inline int bmm_size_cmp(size_t const n, size_t const k) {
  return n < k ? -1 : n > k ? 1 : 0;
}

/// The call `bmm_size_min(n, k)` returns the lesser of `n` and `k`.
/// This is analogous to `fmin`.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_min(size_t const n, size_t const k) {
  return n < k ? n : k;
}

/// The call `bmm_size_max(n, k)` returns the greater of `n` and `k`.
/// This is analogous to `fmax`.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_max(size_t const n, size_t const k) {
  return n > k ? n : k;
}

/// The call `bmm_size_pow(n, k)` returns `n` raised to the power `k`.
/// This is analogous to `pow`.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_pow(size_t const n, size_t const k) {
  size_t m = 1;

  for (size_t i = 0; i < k; ++i)
    m *= n;

  return m;
}

/// The call `bmm_size_identity(n)` returns `n`.
/// This is analogous to `bmm_fp_identity`.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_identity(size_t const n) {
  return n;
}

/// The call `bmm_size_constant(n, k)` returns `n`.
/// This is analogous to `bmm_fp_constant`.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_constant(size_t const n,
    __attribute__ ((__unused__)) size_t const k) {
  return n;
}

/// The call `bmm_size_zero(n)` returns `0`.
/// This is analogous to `bmm_fp_zero`.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_zero(__attribute__ ((__unused__)) size_t const n) {
  return 0;
}

/// The call `bmm_size_one(n)` returns `1`.
/// This is analogous to `bmm_fp_one`.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_one(__attribute__ ((__unused__)) size_t const n) {
  return 1;
}

/// The call `bmm_size_midpoint(n, k)`
/// returns the arithmetic mean of `n` and `k`.
/// This is analogous to `bmm_fp_midpoint`.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_midpoint(size_t const n, size_t const k) {
  // Note that `(n + k) / 2` could wrap and
  // `n < k ? n + (k - n) / 2 : k + (n - k) / 2` could have bad performance.
  return n / 2 + k / 2 + (n % 2 + k % 2) / 2;
}

/// The call `bmm_size_sq(n)` returns `n` squared.
/// This is analogous to `bmm_fp_sq`.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_sq(size_t const n) {
  return n * n;
}

/// The call `bmm_size_firt(n, k)` returns the floor of the `k`th root of `n`.
/// This is analogous to `bmm_fp_rt`.
/// Note that the result may be wrong for large arguments.
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

/// The call `bmm_size_cirt(n, k)` returns the ceiling of the `k`th root of `n`.
/// This is analogous to `bmm_fp_rt`.
/// Note that the result may be wrong for large arguments.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_cirt(size_t const n, size_t const k) {
  return n <= 1 ? n : bmm_size_firt(n - 1, k) + 1;
}

/// The call `bmm_size_flog(n, k)`
/// returns the floor of the base `k` logarithm of `n`.
/// This is analogous to `bmm_fp_log`.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_flog(size_t n, size_t const k) {
  // These do not work because of the attributes.
  // dynamic_assert(n <= 0, "invalid argument");
  // dynamic_assert(k <= 1, "invalid base");

  size_t m = 0;

  while (n >= k) {
    n /= k;
    ++m;
  }

  return m;
}

/// The call `bmm_size_clog(n, k)`
/// returns the ceiling of the base `k` logarithm of `n`.
/// This is analogous to `bmm_fp_log`.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_clog(size_t const n, size_t const k) {
  // These do not work because of the attributes.
  // dynamic_assert(n <= 0, "invalid argument");
  // dynamic_assert(k <= 1, "invalid base");

  return n <= 1 ? 0 : bmm_size_flog(n - 1, k) + 1;
}

/// The call `bmm_size_uclamp(n, b)` returns
///
/// * `n` if `0 <= n < b` and
/// * `b - 1` if `n >= b`.
///
/// This is analogous to `bmm_fp_uclamp`.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_uclamp(size_t const n, size_t const b) {
  return n >= b ? b - 1 : n;
}

/// The call `z = bmm_size_wrap(n, a, b)`
/// solves the periodic equation `z == n - a + k * a` for `z`,
/// where `a <= z < b` and `k` is some integer.
/// This is analogous to `bmm_fp_wrap`.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_wrap(size_t const n, size_t const a, size_t const b) {
  size_t const c = b - a;

  // TODO This is wrong.
  return (n + c - a) % c + a;
}

/// The call `z = bmm_fp_uwrap(n, b)`
/// solves the periodic equation `z == n + k * b` for `z`,
/// where `0 <= z < b` and `k` is some integer.
/// This is analogous to `bmm_fp_uwrap`.
/// The `u` prefix means unsigned or unsymmetric (asymmetric).
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_uwrap(size_t const n, size_t const b) {
  return n % b;
}

/// The call `bmm_size_uinc(n, b)`
/// is equivalent to `bmm_size_uwrap(n + 1, b)` without wrapping.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_uinc(size_t const n, size_t const b) {
  return (n + 1) % b;
}

/// The call `bmm_size_inc(n, a, b)`
/// is equivalent to `bmm_size_wrap(n + 1, a, b)` without wrapping.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_inc(size_t const n, size_t const a, size_t const b) {
  size_t const c = b - a;

  // TODO This is wrong.
  return (n + c - a + 1) % c + a;
}

/// The call `bmm_size_udec(n, b)`
/// is equivalent to `bmm_size_uwrap(n - 1, b)` without wrapping.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_udec(size_t const n, size_t const b) {
  return (n + b - 1) % b;
}

/// The call `bmm_size_dec(n, a, b)`
/// is equivalent to `bmm_size_wrap(n - 1, a, b)` without wrapping.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_dec(size_t const n, size_t const a, size_t const b) {
  size_t const c = b - a;

  // TODO This is wrong.
  return (n + c - a - 1) % c + a;
}

/// The call `bmm_size_sum(n, k)`
/// returns the sum of the array `n` of length `k`.
__attribute__ ((__pure__))
inline size_t bmm_size_sum(size_t const* const n, size_t const k) {
  size_t m = 0;

  for (size_t i = 0; i < k; ++i)
    m += n[i];

  return m;
}

/// The call `bmm_size_prod(n, k)`
/// returns the product of the array `n` of length `k`.
__attribute__ ((__pure__))
inline size_t bmm_size_prod(size_t const* const n, size_t const k) {
  size_t m = 1;

  for (size_t i = 0; i < k; ++i)
    m *= n[i];

  return m;
}

/// The call `bmm_size_hc(pij, i, ndim, nper)`
/// sets the index vector `pij` to the index `i`
/// in a hypercube with dimension `ndim` and side length `nper`.
__attribute__ ((__nonnull__))
inline void bmm_size_hc(size_t* const pij,
    size_t i, size_t const ndim, size_t const nper) {
  for (size_t idim = 0; idim < ndim; ++idim) {
    bmm_size_div_t const qr = bmm_size_div(i, nper);

    i = qr.quot;
    pij[ndim - 1 - idim] = qr.rem;
  }
}

/// The call `bmm_size_unhc(ij, ndim, nper)`
/// returns the index of the index vector `ij`
/// in a hypercube with dimension `ndim` and side length `nper`.
__attribute__ ((__nonnull__, __pure__))
inline size_t bmm_size_unhc(size_t const* const ij,
    size_t const ndim, size_t const nper) {
  size_t i = 0;

  for (size_t idim = 0; idim < ndim; ++idim) {
    i *= nper;
    i += ij[idim];
  }

  return i;
}

/// The call `bmm_size_hc(pij, i, ndim, nper)`
/// sets the index vector `pij` to the index `i`
/// in a hypercuboid with dimension `ndim` and side lengths `nper`.
__attribute__ ((__nonnull__))
inline void bmm_size_hcd(size_t* restrict const pij,
    size_t i, size_t const ndim, size_t const* restrict const nper) {
  for (size_t idim = 0; idim < ndim; ++idim) {
    bmm_size_div_t const qr = bmm_size_div(i, nper[ndim - 1 - idim]);

    i = qr.quot;
    pij[ndim - 1 - idim] = qr.rem;
  }
}

/// The call `bmm_size_unhc(ij, ndim, nper)`
/// returns the index of the index vector `ij`
/// in a hypercuboid with dimension `ndim` and side lengths `nper`.
__attribute__ ((__nonnull__, __pure__))
inline size_t bmm_size_unhcd(size_t const* restrict const ij,
    size_t const ndim, size_t const* restrict const nper) {
  size_t i = 0;

  for (size_t idim = 0; idim < ndim; ++idim) {
    i *= nper[idim];
    i += ij[idim];
  }

  return i;
}

#endif
