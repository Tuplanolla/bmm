/// Common operations for ordered types.

#include <stdbool.h>
#include <stddef.h>

#include "cpp.h"
#include "ext.h"

/// The call `bmm_cmp(x, y)`
/// returns the comparison of `x` and `y`, which is
///
/// * `-1` if `x < y`,
/// * `1` if `x > y` and
/// * `0` otherwise.
///
/// This is useful with `bmm_hsort` for example.
__attribute__ ((__const__, __pure__))
inline int type(bmm_cmp, A)(A const x, A const y) {
  return x < y ? -1 : x > y ? 1 : 0;
}

/// The call `bmm_min(x, y)`
/// returns the lesser of `x` and `y`.
/// If `x` or `y` are not numbers, the behavior is undefined.
__attribute__ ((__const__, __pure__))
inline A type(bmm_min, A)(A const x, A const y) {
  return x < y ? x : y;
}

/// The call `bmm_max(x, y)`
/// returns the lesser of `x` and `y`.
/// If `x` or `y` are not numbers, the behavior is undefined.
__attribute__ ((__const__, __pure__))
inline A type(bmm_max, A)(A const x, A const y) {
  return x > y ? x : y;
}

// TODO Ha ha! Almost everything!

#ifndef EVERYTHING
#define EVERYTHING

/// The call `bmm_size_pow(x, e)`
/// returns `x` to the power of `e`.
/// This is analogous to `bmm_fp_pow`.
/// Overflows are impossible internally but possible externally.
__attribute__ ((__const__, __deprecated__, __pure__))
inline size_t bmm_size_pow(size_t const x, size_t const e) {
  size_t y = 1;

  size_t m = x;
  for (size_t p = e; p > 0; m *= m, p >>= 1)
    if ((p & 1) != 0)
      y *= m;

  return y;

  // The following implementation is less complicated,
  // but slower for large powers.
  // size_t y = 1;
  //
  // for (size_t i = 0; i < e; ++i)
  //   y *= x;
  //
  // return y;
}

/// The call `bmm_size_min(n, k)` returns the lesser of `n` and `k`.
/// This is analogous to `fmin`.
/// Overflows are impossible both internally and externally.
__attribute__ ((__const__, __deprecated__, __pure__))
inline size_t bmm_size_min(size_t const n, size_t const k) {
  return BMM_MIN(n, k);
}

/// The call `bmm_size_max(n, k)` returns the greater of `n` and `k`.
/// This is analogous to `fmax`.
/// Overflows are impossible both internally and externally.
__attribute__ ((__const__, __deprecated__, __pure__))
inline size_t bmm_size_max(size_t const n, size_t const k) {
  return BMM_MAX(n, k);
}

/// The call `bmm_size_identity(n)` returns `n`.
/// This is analogous to `bmm_fp_identity`.
/// Overflows are impossible both internally and externally.
__attribute__ ((__const__, __deprecated__, __pure__))
inline size_t bmm_size_identity(size_t const n) {
  return n;
}

/// The call `bmm_size_constant(n, k)` returns `n`.
/// This is analogous to `bmm_fp_constant`.
/// Overflows are impossible both internally and externally.
__attribute__ ((__const__, __deprecated__, __pure__))
inline size_t bmm_size_constant(size_t const n,
    __attribute__ ((__unused__)) size_t const k) {
  return n;
}

/// The call `bmm_size_zero(n)` returns `0`.
/// This is analogous to `bmm_fp_zero`.
/// Overflows are impossible both internally and externally.
__attribute__ ((__const__, __deprecated__, __pure__))
inline size_t bmm_size_zero(__attribute__ ((__unused__)) size_t const n) {
  return 0;
}

/// The call `bmm_size_one(n)` returns `1`.
/// This is analogous to `bmm_fp_one`.
/// Overflows are impossible both internally and externally.
__attribute__ ((__const__, __deprecated__, __pure__))
inline size_t bmm_size_one(__attribute__ ((__unused__)) size_t const n) {
  return 1;
}

/// The call `bmm_size_midpoint(n, k)`
/// returns the arithmetic mean of `n` and `k`.
/// This is analogous to `bmm_fp_midpoint`.
/// Overflows are impossible both internally and externally.
__attribute__ ((__const__, __deprecated__, __pure__))
inline size_t bmm_size_midpoint(size_t const n, size_t const k) {
  return n / 2 + k / 2 + (n % 2 + k % 2) / 2;

  // The following implementation requires fewer operations,
  // but is slower for branch prediction.
  // return n < k ? n + (k - n) / 2 : k + (n - k) / 2;

  // The following implementation is less complicated,
  // but susceptible to overflowing.
  // return (n + k) / 2;
}

/// The call `bmm_size_sq(n)` returns `n` squared.
/// This is analogous to `bmm_fp_sq`.
/// Overflows are impossible internally but possible externally.
__attribute__ ((__const__, __deprecated__, __pure__))
inline size_t bmm_size_sq(size_t const n) {
  return n * n;
}

/// The call `bmm_size_cb(n)` returns `n` cubed.
/// This is analogous to `bmm_fp_cb`.
/// Overflows are impossible internally but possible externally.
__attribute__ ((__const__, __deprecated__, __pure__))
inline size_t bmm_size_cb(size_t const n) {
  return n * n * n;
}

/// The call `bmm_size_flog(n, k)`
/// returns the floor of the base `k` logarithm of `n`.
/// If `n == 0` or `k < 1`, the behavior is undefined.
/// This is analogous to `bmm_fp_log`.
/// Overflows are impossible both internally and externally.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
__attribute__ ((__deprecated__))
inline size_t bmm_size_flog(size_t n, size_t const k) {
  dynamic_assert(n > 0, "Invalid argument");
  dynamic_assert(k > 1, "Invalid base");

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
/// Overflows are impossible both internally and externally.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
__attribute__ ((__deprecated__))
inline size_t bmm_size_clog(size_t const n, size_t const k) {
  dynamic_assert(n > 0, "Invalid argument");
  dynamic_assert(k > 1, "Invalid base");

  return n <= 1 ? 0 : bmm_size_flog(n - 1, k) + 1;
}

/// The call `bmm_size_firt(n, k)`
/// returns the floor of the `k`th root of `n`.
/// This is analogous to `bmm_fp_rt`.
/// Overflows are possible internally even though they should not be.
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
/// Overflows are possible internally even though they should not be.
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
/// If `b <= 0`, the behavior is undefined.
/// Overflows are impossible both internally and externally.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline size_t bmm_size_uclamp(size_t const n, size_t const b) {
  dynamic_assert(b > 0, "Invalid argument");

  return n >= b ? b - 1 : n;
}

/// The call `bmm_size_uinc(n, b)`
/// is equivalent to `bmm_size_inc(n, 0, b)`.
/// If `b <= 0`, the behavior is undefined.
/// Overflows are possible both internally and externally.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline size_t bmm_size_uinc(size_t const n, size_t const b) {
  dynamic_assert(b > 0, "Invalid argument");

  return (n + 1) % b;
}

/// The call `bmm_size_inc(n, a, b)`
/// is equivalent to `type(bmm_wrap, size_t)(n + 1, a, b)` without wrapping.
/// If `b <= a`, the behavior is undefined.
/// Overflows are possible both internally and externally.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline size_t bmm_size_inc(size_t const n, size_t const a, size_t const b) {
  dynamic_assert(b > a, "Invalid argument");

  size_t const c = b - a;

  return (n % c + c - a % c + 1) % c + a;
}

/// The call `bmm_size_udec(n, b)`
/// is equivalent to `bmm_size_dec(n, 0, b)`.
/// If `b <= 0`, the behavior is undefined.
/// Overflows are possible both internally and externally.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline size_t bmm_size_udec(size_t const n, size_t const b) {
  dynamic_assert(b > 0, "Invalid argument");

  return (n + b - 1) % b;
}

/// The call `bmm_size_dec(n, a, b)`
/// is equivalent to `type(bmm_wrap, size_t)(n - 1, a, b)` without wrapping.
/// If `b <= a`, the behavior is undefined.
/// Overflows are possible both internally and externally.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline size_t bmm_size_dec(size_t const n, size_t const a, size_t const b) {
  dynamic_assert(b > a, "Invalid argument");

  size_t const c = b - a;

  return (n % c + c - a % c - 1) % c + a;
}

/// The call `bmm_size_fact(n, k)`
/// returns the `k`-factorial of `n`.
/// Overflows are impossible internally but possible externally.
__attribute__ ((__const__, __deprecated__, __pure__))
inline size_t bmm_size_fact(size_t const n, size_t const k) {
  size_t m = 1;

  for (size_t i = n; i > 1; i -= k)
    m *= i;

  return m;
}

/// The call `bmm_size_tri(n)`
/// returns the `n`th triangular number.
/// Overflows are impossible internally but possible externally.
__attribute__ ((__const__, __pure__))
inline size_t bmm_size_tri(size_t const n) {
  // return bmm_size_choose(n + 1, 2);
  return n * (n + 1) / 2;
}

/// The call `bmm_size_sum(n, k)`
/// returns the sum of the array `n` of length `k`.
/// Overflows are impossible internally but possible externally.
__attribute__ ((__pure__))
inline size_t bmm_size_sum(size_t const *const n, size_t const k) {
  size_t m = 0;

  for (size_t i = 0; i < k; ++i)
    m += n[i];

  return m;
}

/// The call `bmm_size_prod(n, k)`
/// returns the product of the array `n` of length `k`.
/// Overflows are impossible internally but possible externally.
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
__attribute__ ((__deprecated__, __nonnull__ (1, 2)))
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
__attribute__ ((__deprecated__, __nonnull__ (1, 2)))
inline size_t bmm_size_rfold(size_t (*const f)(size_t, size_t, void *),
    size_t const *restrict const x, size_t const n,
    size_t z, void *restrict const ptr) {
  for (size_t i = 0; i < n; ++i)
    z = f(x[n - 1 - i], z, ptr);

  return z;
}

#endif
