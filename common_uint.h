/// Common operations for unsigned integer types.

#include <stddef.h>

#include "ext.h"

/// The call `bmm_quotrem(x, y)`
/// returns the quotient and remainder of `x` divided by `y`
/// in `qr` such that `qr.quot * y + qr.rem == x` and `qr.rem >= 0`.
/// If `y == 0`, the behavior is undefined.
/// Overflows are impossible both internally and externally.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline type(bmm_quotrem_t, A) type(bmm_quotrem, A)(A const x, A const y) {
  dynamic_assert(y != 0, "Invalid argument");

  return (type(bmm_quotrem_t, A)) {.quot = x / y, .rem = x % y};
}

/// The call `bmm_abs(x)`
/// returns the absolute value of `x`.
/// Overflows are impossible both internally and externally.
__attribute__ ((__const__, __pure__))
inline A type(bmm_abs, A)(A const x) {
  return x;
}

/// The call `bmm_wrap(x, a, b)`
/// finds such `y` that `a <= y < b`
/// by shifting `x` by the appropriate number of multiples of `b - a`.
/// If `b <= a`, the behavior is undefined.
/// Overflows are impossible both internally and externally.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A type(bmm_wrap, A)(A const x, A const a, A const b) {
  dynamic_assert(b > a, "Invalid argument");

  A const c = b - a;
  A const r = x % c;
  A const s = a % c;

  return (r >= s ? r - s : c - (s - r)) + a;

  // The following implementation is easier to understand,
  // but susceptible to overflowing.
  // A const c = b - a;
  //
  // return (x - a) % c + a;

  // The following implementation is easier to understand,
  // but slower (linear instead of constant).
  // A const c = b - a;
  //
  // A y = x;
  //
  // if (y < a)
  //   do
  //     y += c;
  //   while (y < a);
  // else if (y >= b)
  //   do
  //     y -= c;
  //   while (y >= b);
  //
  // return y;
}

/// The call `bmm_uwrap(x, b)`
/// is equivalent to `bmm_wrap(x, 0, b)`.
/// If `b <= 0`, the behavior is undefined.
/// Overflows are impossible both internally and externally.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A type(bmm_uwrap, A)(A const x, A const b) {
  dynamic_assert(b > 0, "Invalid argument");

  return x % b;
}

/// The call `bmm_uclamp(x, b)`
/// finds such `y` that `0 <= y <= b`
/// by shifting `x` by the smallest possible amount.
/// Overflows are impossible both internally and externally.
__attribute__ ((__const__, __pure__))
inline A type(bmm_uclamp, A)(A const x, A const b) {
  return x > b ? b : x;
}

// TODO Do we even need these? How about `bmm_wrapadd(x, y, a, b)`?

/// The call `bmm_uinc(x, b)`
/// is equivalent to `bmm_inc(x, 0, b)`.
/// If `b <= 0`, the behavior is undefined.
/// Overflows are possible both internally and externally.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A type(bmm_uinc, A)(A const x, A const b) {
  dynamic_assert(b > 0, "Invalid argument");

  return (x + 1) % b;
}

/// The call `bmm_inc(x, a, b)`
/// is equivalent to `type(bmm_wrap, A)(x + 1, a, b)` without wrapping.
/// If `b <= a`, the behavior is undefined.
/// Overflows are possible both internally and externally.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A type(bmm_inc, A)(A const x, A const a, A const b) {
  dynamic_assert(b > a, "Invalid argument");

  A const c = b - a;

  return (x % c + c - a % c + 1) % c + a;
}

/// The call `bmm_udec(x, b)`
/// is equivalent to `bmm_dec(x, 0, b)`.
/// If `b <= 0`, the behavior is undefined.
/// Overflows are possible both internally and externally.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A type(bmm_udec, A)(A const x, A const b) {
  dynamic_assert(b > 0, "Invalid argument");

  return (x + b - 1) % b;
}

/// The call `bmm_dec(x, a, b)`
/// is equivalent to `type(bmm_wrap, A)(x - 1, a, b)` without wrapping.
/// If `b <= a`, the behavior is undefined.
/// Overflows are possible both internally and externally.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A type(bmm_dec, A)(A const x, A const a, A const b) {
  dynamic_assert(b > a, "Invalid argument");

  A const c = b - a;

  return (x % c + c - a % c - 1) % c + a;
}

/// The call `bmm_fact(x)`
/// returns the factorial of `x`.
/// Overflows are impossible internally but possible externally.
__attribute__ ((__const__, __pure__))
inline A type(bmm_fact, A)(A const x) {
  if (x <= 1)
    return 1;

  A y = x;

  A z = y;
  while (z > 1) {
    --z;
    y *= z;
  }

  return y;
}

/// The call `bmm_multfact(x, m)`
/// returns the multifactorial of `x` with the multiplicity `m`.
/// If `m == 0`, the behavior is undefined.
/// Overflows are impossible internally but possible externally.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A type(bmm_multfact, A)(A const x, A const m) {
  dynamic_assert(m > 0, "Invalid argument");

  if (x <= 1)
    return 1;

  A y = x;

  A z = y;
  while (z > m) {
    z -= m;
    y *= z;
  }

  return y;
}

// TODO Check and test this.

/// The call `bmm_tri(x)`
/// returns the `x`th triangular number.
/// Overflows are possible internally even though they should not be.
__attribute__ ((__const__, __pure__))
inline A type(bmm_tri, A)(A const x) {
  // return type(bmm_choose, A)(x + 1, 2);
  return x * (x + 1) / 2;
}

/// The call `bmm_tamean2(x, y)`
/// returns the truncated arithmetic mean of `x` and `y`.
/// Overflows are impossible both internally and externally.
__attribute__ ((__const__, __pure__))
inline A type(bmm_tamean2, A)(A const x, A const y) {
  return (x / 2 + y / 2) + (x % 2) * (y % 2);
}

/// The call `bmm_famean2(x, y)`
/// returns the floored arithmetic mean of `x` and `y`.
/// Overflows are impossible both internally and externally.
__attribute__ ((__const__, __pure__))
inline A type(bmm_famean2, A)(A const x, A const y) {
  return (x / 2 + y / 2) + (x % 2) * (y % 2);
}

/// The call `bmm_flog(x, b)`
/// returns the floored base `b` logarithm of `x`.
/// If `x == 0` or `b <= 1`, the behavior is undefined.
/// Overflows are impossible both internally and externally.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A type(bmm_flog, A)(A const x, A const b) {
  dynamic_assert(x > 0, "Invalid argument");
  dynamic_assert(b > 1, "Invalid argument");

  A y = 0;

  A z = x;
  while (z >= b) {
    z /= b;
    ++y;
  }

  return y;
}

/// The call `bmm_clog(x, b)`
/// returns the ceilinged base `b` logarithm of `x`.
/// If `x == 0` or `b <= 1`, the behavior is undefined.
/// Overflows are impossible both internally and externally.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A type(bmm_clog, A)(A const x, A const b) {
  dynamic_assert(x > 0, "Invalid argument");
  dynamic_assert(b > 1, "Invalid argument");

  return x == 1 ? 0 : type(bmm_flog, A)(x - 1, b) + 1;
}

// TODO Check and test these.

/// The call `bmm_firt(n, k)`
/// returns the floored root `k` of `n`.
/// Overflows are possible internally even though they should not be.
__attribute__ ((__const__, __pure__))
inline A type(bmm_firt, A)(A const n, A const k) {
  if (n <= 1)
    return n;
  else {
    A const p = k - 1;
    A r = n + 1;
    A m = n;

    while (m < r) {
      r = m;
      m = (p * r + n / type(bmm_power, A)(r, p)) / k;
    }

    return r;
  }
}

/// The call `bmm_cirt(n, k)`
/// returns the ceilinged root `k` of `n`.
/// Overflows are possible internally even though they should not be.
__attribute__ ((__const__, __pure__))
inline A type(bmm_cirt, A)(A const n, A const k) {
  return n <= 1 ? n : type(bmm_firt, A)(n - 1, k) + 1;
}

/// The call `bmm_hc(pij, i, ndim, nper)`
/// sets the index vector `pij` to the index `i`
/// in a hypercube with dimension `ndim` and side length `nper`.
/// Overflows are impossible both internally and externally.
__attribute__ ((__nonnull__))
inline void type(bmm_hc, A)(A *const pij, A const i,
    size_t const ndim, A const nper) {
  type(bmm_quotrem_t, A) dm = {.quot = i};
  for (size_t idim = 0; idim < ndim; ++idim) {
    dm = type(bmm_quotrem, A)(dm.quot, nper);

    pij[ndim - 1 - idim] = dm.rem;
  }
}

/// The call `bmm_unhc(ij, ndim, nper)`
/// returns the index of the index vector `ij`
/// in a hypercube with dimension `ndim` and side length `nper`.
/// Overflows are impossible internally but possible externally.
__attribute__ ((__nonnull__, __pure__))
inline A type(bmm_unhc, A)(A const *const ij,
    size_t const ndim, A const nper) {
  A i = 0;

  for (size_t idim = 0; idim < ndim; ++idim) {
    i *= nper;
    i += ij[idim];
  }

  return i;
}

/// The call `bmm_hc(pij, i, ndim, nper)`
/// sets the index vector `pij` to the index `i`
/// in a hypercuboid with dimension `ndim` and side lengths `nper`.
/// Overflows are impossible both internally and externally.
__attribute__ ((__nonnull__))
inline void type(bmm_hcd, A)(A *restrict const pij, A const i,
    size_t const ndim, A const *restrict const nper) {
  type(bmm_quotrem_t, A) dm = {.quot = i};
  for (size_t idim = 0; idim < ndim; ++idim) {
    dm = type(bmm_quotrem, A)(dm.quot, nper[ndim - 1 - idim]);

    pij[ndim - 1 - idim] = dm.rem;
  }

  // The following implementation is suitable for loop fusion,
  // but less reliable.
  // size_t *const buf = alloca(ndim * sizeof *buf);
  //
  // type(bmm_quotrem_t, A) dm = {.quot = i};
  // for (size_t idim = 0; idim < ndim; ++idim) {
  //   dm = type(bmm_quotrem, A)(dm.quot, nper[ndim - 1 - idim]);
  //
  //   buf[ndim - 1 - idim] = dm.rem;
  // }
  //
  // for (size_t idim = 0; idim < ndim; ++idim)
  //   pij[idim] = buf[idim];

  // The following implementation is suitable for loop fusion,
  // but slower (quadratic instead of linear).
  // for (size_t idim = 0; idim < ndim; ++idim) {
  //   type(bmm_quotrem_t, A) dm = {.quot = i};
  //   for (size_t jdim = 0; jdim < ndim - idim; ++jdim)
  //     dm = type(bmm_quotrem, A)(dm.quot, nper[ndim - 1 - jdim]);
  //
  //   pij[idim] = dm.rem;
  // }
}

/// The call `bmm_unhc(ij, ndim, nper)`
/// returns the index of the index vector `ij`
/// in a hypercuboid with dimension `ndim` and side lengths `nper`.
/// Overflows are impossible internally but possible externally.
__attribute__ ((__nonnull__, __pure__))
inline A type(bmm_unhcd, A)(A const *restrict const ij,
    size_t const ndim, A const *restrict const nper) {
  A i = 0;

  for (size_t idim = 0; idim < ndim; ++idim) {
    i *= nper[idim];
    i += ij[idim];
  }

  return i;
}
