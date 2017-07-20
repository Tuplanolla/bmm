#include <stdbool.h>
#include <stddef.h>

#include "cpp.h"
#include "ext.h"

// TODO Ha ha! Everything!

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
/// Overflows are handled appropriately.
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

/// The call `bmm_size_swap(x, y)`
/// exchanges `x` with `y`.
__attribute__ ((__nonnull__))
inline void bmm_size_swap(size_t *restrict const x, size_t *restrict const y) {
  size_t const tmp = *x;
  *x = *y;
  *y = tmp;
}

/// The call `bmm_size_even(n)`
/// checks whether `n` is even.
__attribute__ ((__const__, __pure__))
inline bool bmm_size_even(size_t const n) {
  return n % 2 == 0;
}

/// The call `bmm_size_odd(n)`
/// checks whether `n` is odd.
__attribute__ ((__const__, __pure__))
inline bool bmm_size_odd(size_t const n) {
  return n % 2 == 1;
}

/// The call `bmm_size_min(n, k)` returns the lesser of `n` and `k`.
/// This is analogous to `fmin`.
/// Overflows are handled appropriately.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_min(size_t const n, size_t const k) {
  return n < k ? n : k;
}

/// The call `bmm_size_max(n, k)` returns the greater of `n` and `k`.
/// This is analogous to `fmax`.
/// Overflows are handled appropriately.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_max(size_t const n, size_t const k) {
  return n > k ? n : k;
}

/// The call `bmm_size_identity(n)` returns `n`.
/// This is analogous to `bmm_fp_identity`.
/// Overflows are handled appropriately.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_identity(size_t const n) {
  return n;
}

/// The call `bmm_size_constant(n, k)` returns `n`.
/// This is analogous to `bmm_fp_constant`.
/// Overflows are handled appropriately.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_constant(size_t const n,
    __attribute__ ((__unused__)) size_t const k) {
  return n;
}

/// The call `bmm_size_zero(n)` returns `0`.
/// This is analogous to `bmm_fp_zero`.
/// Overflows are handled appropriately.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_zero(__attribute__ ((__unused__)) size_t const n) {
  return 0;
}

/// The call `bmm_size_one(n)` returns `1`.
/// This is analogous to `bmm_fp_one`.
/// Overflows are handled appropriately.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_one(__attribute__ ((__unused__)) size_t const n) {
  return 1;
}

/// The call `bmm_size_midpoint(n, k)`
/// returns the arithmetic mean of `n` and `k`.
/// This is analogous to `bmm_fp_midpoint`.
/// Overflows are handled appropriately.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_midpoint(size_t const n, size_t const k) {
  return n / 2 + k / 2 + (n % 2 + k % 2) / 2;

  // The following implementation is less laborious,
  // but slower for branch prediction.
  // return n < k ? n + (k - n) / 2 : k + (n - k) / 2;

  // The following implementation is less complicated,
  // but susceptible to overflowing.
  // return (n + k) / 2;
}

/// The call `bmm_size_sq(n)` returns `n` squared.
/// This is analogous to `bmm_fp_sq`.
/// Overflows are not handled appropriately.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_sq(size_t const n) {
  return n * n;
}

/// The call `bmm_size_cb(n)` returns `n` cubed.
/// This is analogous to `bmm_fp_cb`.
/// Overflows are not handled appropriately.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_cb(size_t const n) {
  return n * n * n;
}

/// The call `bmm_size_flog(n, k)`
/// returns the floor of the base `k` logarithm of `n`.
/// If `n == 0` or `k < 1`, the behavior is undefined.
/// This is analogous to `bmm_fp_log`.
/// Overflows are handled appropriately.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline size_t bmm_size_flog(size_t n, size_t const k) {
#ifdef DEBUG
  // These do not work together with the attributes.
  dynamic_assert(n > 0, "Invalid argument");
  dynamic_assert(k > 1, "Invalid base");
#endif

  size_t m = 0;

  while (n >= k) {
    n /= k;
    ++m;
  }

  return m;
}

/// The call `bmm_size_clog(n, k)`
/// returns the ceiling of the base `k` logarithm of `n`.
/// If `n == 0` or `k < 1`, the behavior is undefined.
/// This is analogous to `bmm_fp_log`.
/// Overflows are handled appropriately.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline size_t bmm_size_clog(size_t const n, size_t const k) {
#ifdef DEBUG
  // These do not work together with the attributes.
  dynamic_assert(n > 0, "Invalid argument");
  dynamic_assert(k > 1, "Invalid base");
#endif

  return n <= 1 ? 0 : bmm_size_flog(n - 1, k) + 1;
}

/// The call `bmm_size_pow(x, e)`
/// returns `x` to the power of `e`.
/// This is analogous to `bmm_fp_pow`.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_pow(size_t const x, size_t const e) {
  size_t y = 1;

  size_t m = x;
  for (size_t p = e; p > 0; m *= m, p >>= 1)
    if ((p & 1) != 0)
      y *= m;

  return y;

  // The following implementation is less complicated,
  // but slower for large powers.
  // size_t m = 1;
  //
  // for (size_t i = 0; i < k; ++i)
  //   m *= n;
  //
  // return m;
}

/// The call `bmm_size_firt(n, k)`
/// returns the floor of the `k`th root of `n`.
/// This is analogous to `bmm_fp_rt`.
/// Overflows are not handled appropriately.
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

/// The call `bmm_size_cirt(n, k)`
/// returns the ceiling of the `k`th root of `n`.
/// This is analogous to `bmm_fp_rt`.
/// Overflows are not handled appropriately.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_cirt(size_t const n, size_t const k) {
  return n <= 1 ? n : bmm_size_firt(n - 1, k) + 1;
}

/// The call `bmm_size_uclamp(n, b)` returns
///
/// * `n` if `0 <= n < b` and
/// * `b - 1` if `n >= b`.
///
/// This is analogous to `bmm_fp_uclamp`.
/// Overflows are handled appropriately.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_uclamp(size_t const n, size_t const b) {
  return n >= b ? b - 1 : n;
}

/// The call `m = bmm_size_wrap(n, a, b)`
/// solves the periodic equation `m == n - a + k * a` for `m`,
/// where `a <= m < b` and `k` is some integer.
/// This is analogous to `bmm_fp_wrap`.
/// Overflows are handled appropriately.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_wrap(size_t const n, size_t const a, size_t const b) {
  size_t const c = b - a;

  return (n % c + c - a % c) % c + a;

  // The following implementation is a lot slower, but easier to understand.
  // size_t const c = b - a;
  //
  // size_t k = n;
  //
  // if (k < a)
  //   do
  //     k += c;
  //   while (k < a);
  // else if (k >= b)
  //   do
  //     k -= c;
  //   while (k >= b);
  //
  // return k;
}

/// The call `m = bmm_size_uwrap(n, b)`
/// solves the periodic equation `m == n + k * b` for `m`,
/// where `0 <= m < b` and `k` is some integer.
/// This is analogous to `bmm_fp_uwrap`.
/// The `u` prefix means unsigned or unsymmetric (asymmetric).
/// Overflows are handled appropriately.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_uwrap(size_t const n, size_t const b) {
  return n % b;
}

/// The call `bmm_size_uinc(n, b)`
/// is equivalent to `bmm_size_inc(n, 0, b)`.
/// Overflows are handled appropriately.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_uinc(size_t const n, size_t const b) {
  return (n + 1) % b;
}

/// The call `bmm_size_inc(n, a, b)`
/// is equivalent to `bmm_size_wrap(n + 1, a, b)` without wrapping.
/// Overflows are handled appropriately.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_inc(size_t const n, size_t const a, size_t const b) {
  size_t const c = b - a;

  return (n % c + c - a % c + 1) % c + a;
}

/// The call `bmm_size_udec(n, b)`
/// is equivalent to `bmm_size_dec(n, 0, b)`.
/// Overflows are handled appropriately.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_udec(size_t const n, size_t const b) {
  return (n + b - 1) % b;
}

/// The call `bmm_size_dec(n, a, b)`
/// is equivalent to `bmm_size_wrap(n - 1, a, b)` without wrapping.
/// Overflows are handled appropriately.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_dec(size_t const n, size_t const a, size_t const b) {
  size_t const c = b - a;

  return (n % c + c - a % c - 1) % c + a;
}

/// The call `bmm_size_fact(n, k)`
/// returns the `k`-factorial of `n`.
/// Overflows are not handled appropriately.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_fact(size_t const n, size_t const k) {
  size_t m = 1;

  for (size_t i = n; i > 1; i -= k)
    m *= i;

  return m;
}

/// The call `bmm_size_tri(n)`
/// returns the `n`th triangular number.
/// Overflows are not handled appropriately.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_tri(size_t const n) {
  // return bmm_size_choose(n + 1, 2);
  return n * (n + 1) / 2;
}

/// The call `bmm_size_sum(n, k)`
/// returns the sum of the array `n` of length `k`.
/// Overflows are not handled appropriately.
__attribute__ ((__pure__))
inline size_t bmm_size_sum(size_t const *const n, size_t const k) {
  size_t m = 0;

  for (size_t i = 0; i < k; ++i)
    m += n[i];

  return m;
}

/// The call `bmm_size_prod(n, k)`
/// returns the product of the array `n` of length `k`.
/// Overflows are not handled appropriately.
__attribute__ ((__pure__))
inline size_t bmm_size_prod(size_t const *const n, size_t const k) {
  size_t m = 1;

  for (size_t i = 0; i < k; ++i)
    m *= n[i];

  return m;
}

/// The call `bmm_size_lfold(f, x, n, z, ptr)`
/// folds the procedure `f` over the array `x` of length `n`.
/// by starting from the left with `z`.
__attribute__ ((__nonnull__ (1, 2)))
inline size_t bmm_size_lfold(size_t (*const f)(size_t, size_t, void *),
    size_t const *restrict const x, size_t const n,
    size_t z, void *restrict const ptr) {
  for (size_t i = 0; i < n; ++i)
    z = f(x[i], z, ptr);

  return z;
}

/// The call `bmm_size_rfold(f, x, n, z, ptr)`
/// folds the procedure `f` over the array `x` of length `n`.
/// by starting from the right with `z`.
__attribute__ ((__nonnull__ (1, 2)))
inline size_t bmm_size_rfold(size_t (*const f)(size_t, size_t, void *),
    size_t const *restrict const x, size_t const n,
    size_t z, void *restrict const ptr) {
  for (size_t i = 0; i < n; ++i)
    z = f(x[n - 1 - i], z, ptr);

  return z;
}

/// The call `bmm_size_hc(pij, i, ndim, nper)`
/// sets the index vector `pij` to the index `i`
/// in a hypercube with dimension `ndim` and side length `nper`.
/// Overflows are handled appropriately.
__attribute__ ((__nonnull__))
inline void bmm_size_hc(size_t *const pij,
    size_t const i, size_t const ndim, size_t const nper) {
  bmm_size_div_t qr = {.quot = i, .rem = 0};
  for (size_t idim = 0; idim < ndim; ++idim) {
    qr = bmm_size_div(qr.quot, nper);

    pij[ndim - 1 - idim] = qr.rem;
  }
}

/// The call `bmm_size_unhc(ij, ndim, nper)`
/// returns the index of the index vector `ij`
/// in a hypercube with dimension `ndim` and side length `nper`.
/// Overflows are not handled appropriately.
__attribute__ ((__nonnull__, __pure__))
inline size_t bmm_size_unhc(size_t const *const ij,
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
/// Overflows are handled appropriately.
__attribute__ ((__nonnull__))
inline void bmm_size_hcd(size_t *restrict const pij,
    size_t const i, size_t const ndim, size_t const *restrict const nper) {
  bmm_size_div_t qr = {.quot = i, .rem = 0};
  for (size_t idim = 0; idim < ndim; ++idim) {
    qr = bmm_size_div(qr.quot, nper[ndim - 1 - idim]);

    pij[ndim - 1 - idim] = qr.rem;
  }

  // The following implementation is less reliable,
  // but suitable for loop fusion.
  // size_t *const buf = alloca(ndim * sizeof *buf);
  //
  // bmm_size_div_t qr = {.quot = i, .rem = 0};
  // for (size_t idim = 0; idim < ndim; ++idim) {
  //   qr = bmm_size_div(qr.quot, nper[ndim - 1 - idim]);
  //
  //   buf[ndim - 1 - idim] = qr.rem;
  // }
  //
  // for (size_t idim = 0; idim < ndim; ++idim)
  //   pij[idim] = buf[idim];

  // The following implementation is slower, but suitable for loop fusion.
  // for (size_t idim = 0; idim < ndim; ++idim) {
  //   bmm_size_div_t qr = {.quot = i, .rem = 0};
  //   for (size_t jdim = 0; jdim < ndim - idim; ++jdim)
  //     qr = bmm_size_div(qr.quot, nper[ndim - 1 - jdim]);
  //
  //   pij[idim] = qr.rem;
  // }
}

/// The call `bmm_size_unhc(ij, ndim, nper)`
/// returns the index of the index vector `ij`
/// in a hypercuboid with dimension `ndim` and side lengths `nper`.
/// Overflows are not handled appropriately.
__attribute__ ((__nonnull__, __pure__))
inline size_t bmm_size_unhcd(size_t const *restrict const ij,
    size_t const ndim, size_t const *restrict const nper) {
  size_t i = 0;

  for (size_t idim = 0; idim < ndim; ++idim) {
    i *= nper[idim];
    i += ij[idim];
  }

  return i;
}
