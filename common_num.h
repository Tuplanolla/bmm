/// Common operations for numeric types.

#include <stddef.h>

#include "ext.h"

/// This structure holds the quotient and remainder of division
/// in unspecified order.
typedef struct {
  A quot;
  A rem;
} $(bmm_quotrem_t, A);

// This forward-declaration makes reverse dependencies possible.
inline $(bmm_quotrem_t, A) $(bmm_quotrem, A)(A, A);

/// The call `bmm_quot(x, y)`
/// returns the quotient of `x` divided by `y`.
/// If `y == 0`, the behavior is undefined.
/// Overflows are impossible both internally and externally.
#ifndef DEBUG
__attribute__ ((__const__, __flatten__, __pure__))
#endif
inline A $(bmm_quot, A)(A const x, A const y) {
  return $(bmm_quotrem, A)(x, y).quot;
}

/// The call `bmm_rem(x, y)`
/// returns the remainder of `x` divided by `y`.
/// If `y == 0`, the behavior is undefined.
/// Overflows are impossible both internally and externally.
#ifndef DEBUG
__attribute__ ((__const__, __flatten__, __pure__))
#endif
inline A $(bmm_rem, A)(A const x, A const y) {
  return $(bmm_quotrem, A)(x, y).rem;
}

/// The call `bmm_sgn(x)`
/// returns the sign of `x`.
/// If `x` is not a number, the behavior is undefined.
/// Overflows are impossible both internally and externally.
__attribute__ ((__const__, __flatten__, __pure__))
inline int $(bmm_sgn, A)(A const x) {
  return $(bmm_cmp, A)(x, 0);
}

/// The call `bmm_power(x, e)`
/// returns `x` raised to the power of `e`.
/// If `x == 0` and `e == 0`, the result is one.
/// If `x` is not a number, the behavior is undefined.
/// Overflows are impossible internally but possible externally.
__attribute__ ((__const__, __pure__))
inline A $(bmm_power, A)(A const x, size_t const e) {
  A y = 1;

  A m = x;
  for (size_t p = e; p > 0; m *= m, p >>= 1)
    if ((p & 1) != 0)
      y *= m;

  return y;

  // The following implementation is less complicated,
  // but slower for large powers.
  // A y = 1;
  //
  // for (size_t i = 0; i < e; ++i)
  //   y *= x;
  //
  // return y;
}

/// The call `bmm_sum(x, nmemb)`
/// returns the sum of the array `x` of length `nmemb`.
/// Overflows are impossible internally but possible externally.
__attribute__ ((__pure__))
inline A $(bmm_sum, A)(A const *const x, size_t const nmemb) {
  A y = 0;

  for (size_t i = 0; i < nmemb; ++i)
    y += x[i];

  return y;
}

/// The call `bmm_prod(x, nmemb)`
/// returns the product of the array `x` of length `nmemb`.
/// Overflows are impossible internally but possible externally.
__attribute__ ((__pure__))
inline A $(bmm_prod, A)(A const *const x, size_t const nmemb) {
  A y = 1;

  for (size_t i = 0; i < nmemb; ++i)
    y *= x[i];

  return y;
}
