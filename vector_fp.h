/// The call `zero(ox)`
/// stores into `ox` zero.
/// This is analogous to the constant `0`.
__attribute__ ((__nonnull__))
inline void $(zero, A, D)(A *const ox) {
  for (size_t i = 0; i < D; ++i)
    ox[i] = $(zero, A)();
}

/// The call `add(oz, x, y)`
/// stores into `oz` the sum of `x` and `y`.
/// This is analogous to the binary operator `+`.
__attribute__ ((__nonnull__))
inline void $(add, A, D)(A *restrict const oz,
    A const *restrict const x, A const *restrict const y) {
  for (size_t i = 0; i < D; ++i)
    oz[i] = $(add, A)(x[i], y[i]);
}

/// The call `neg(oy, x)`
/// stores into `oy` the negation of `x`.
/// This is analogous to the unary operator `-`.
__attribute__ ((__nonnull__))
inline void $(neg, A, D)(A *restrict const oy, A const *restrict const x) {
  for (size_t i = 0; i < D; ++i)
    oy[i] = $(neg, A)(x[i]);
}

/// The call `sub(oz, x, y)`
/// stores into `oz` the difference of `x` and `y`.
/// This is analogous to the binary operator `-`.
__attribute__ ((__nonnull__))
inline void $(sub, A, D)(A *restrict const oz,
    A const *restrict const x, A const *restrict const y) {
  for (size_t i = 0; i < D; ++i)
    oz[i] = $(sub, A)(x[i], y[i]);
}

/// The call `one(ox)`
/// stores into `ox` one.
/// This is analogous to the constant `1`.
__attribute__ ((__nonnull__))
inline void $(one, A, D)(A *const ox) {
  for (size_t i = 0; i < D; ++i)
    ox[i] = $(one, A)();
}

/// The call `mul(oz, x, y)`
/// stores into `oz` the product of `x` and `y`.
/// This is analogous to the binary operator `*`.
__attribute__ ((__nonnull__))
inline void $(mul, A, D)(A *restrict const oz,
    A const *restrict const x, A const *restrict const y) {
  for (size_t i = 0; i < D; ++i)
    oz[i] = $(mul, A)(x[i], y[i]);
}

/// The call `recip(oy, x)`
/// stores into `oy` the reciprocal of `x`.
__attribute__ ((__nonnull__))
inline void $(recip, A, D)(A *restrict const oy,
    A const *restrict const x) {
  for (size_t i = 0; i < D; ++i)
    oy[i] = $(recip, A)(x[i]);
}

/// The call `div(oz, x, y)`
/// stores into `oz` the division of `x` and `y`.
/// This is analogous to the binary operator `/`.
__attribute__ ((__nonnull__))
inline void $(div, A, D)(A *restrict const oz,
    A const *restrict const x, A const *restrict const y) {
  for (size_t i = 0; i < D; ++i)
    oz[i] = $(div, A)(x[i], y[i]);
}

/// The call `quott(oz, x, y)`
/// stores into `oz` the truncated quotient of `x` and `y`.
__attribute__ ((__nonnull__))
inline void $(quott, A, D)(A *restrict const oz,
    A const *restrict const x, A const *restrict const y) {
  for (size_t i = 0; i < D; ++i)
    oz[i] = $(quott, A)(x[i], y[i]);
}

/// The call `remt(oz, x, y)`
/// stores into `oz` the truncated remainder of `x` and `y`.
__attribute__ ((__nonnull__))
inline void $(remt, A, D)(A *restrict const oz,
    A const *restrict const x, A const *restrict const y) {
  for (size_t i = 0; i < D; ++i)
    oz[i] = $(remt, A)(x[i], y[i]);
}

/// The call `quote(oz, x, y)`
/// stores into `oz` the Euclidean quotient of `x` and `y`.
__attribute__ ((__nonnull__))
inline void $(quote, A, D)(A *restrict const oz,
    A const *restrict const x, A const *restrict const y) {
  for (size_t i = 0; i < D; ++i)
    oz[i] = $(quote, A)(x[i], y[i]);
}

/// The call `reme(oz, x, y)`
/// stores into `oz` the Euclidean remainder of `x` and `y`.
__attribute__ ((__nonnull__))
inline void $(reme, A, D)(A *restrict const oz,
    A const *restrict const x, A const *restrict const y) {
  for (size_t i = 0; i < D; ++i)
    oz[i] = $(reme, A)(x[i], y[i]);
}

/// The call `add_mut(iox, y)`
/// stores into `iox` the sum of `iox` and `y`.
/// This is analogous to the binary operator `+=`.
__attribute__ ((__nonnull__))
inline void $(add_mut, A, D)(A *restrict const iox,
    A const *restrict const y) {
  for (size_t i = 0; i < D; ++i)
    iox[i] = $(add, A)(iox[i], y[i]);
}

/// The call `neg_mut(iox)`
/// stores into `iox` the negation of `iox`.
__attribute__ ((__nonnull__))
inline void $(neg_mut, A, D)(A *restrict const iox) {
  for (size_t i = 0; i < D; ++i)
    iox[i] = $(neg, A)(iox[i]);
}

/// The call `sub_mut(iox, y)`
/// stores into `iox` the difference of `iox` and `y`.
/// This is analogous to the binary operator `-=`.
__attribute__ ((__nonnull__))
inline void $(sub_mut, A, D)(A *restrict const iox,
    A const *restrict const y) {
  for (size_t i = 0; i < D; ++i)
    iox[i] = $(sub, A)(iox[i], y[i]);
}

/// The call `mul_mut(iox, y)`
/// stores into `iox` the product of `iox` and `y`.
/// This is analogous to the binary operator `*=`.
__attribute__ ((__nonnull__))
inline void $(mul_mut, A, D)(A *restrict const iox,
    A const *restrict const y) {
  for (size_t i = 0; i < D; ++i)
    iox[i] = $(mul, A)(iox[i], y[i]);
}

/// The call `recip_mut(iox)`
/// stores into `iox` the reciprocal of `iox`.
__attribute__ ((__nonnull__))
inline void $(recip_mut, A, D)(A *restrict const iox) {
  for (size_t i = 0; i < D; ++i)
    iox[i] = $(recip, A)(iox[i]);
}

/// The call `div_mut(iox, y)`
/// stores into `iox` the division of `iox` and `y`.
/// This is analogous to the binary operator `/=`.
__attribute__ ((__nonnull__))
inline void $(div_mut, A, D)(A *restrict const iox,
    A const *restrict const y) {
  for (size_t i = 0; i < D; ++i)
    iox[i] = $(div, A)(iox[i], y[i]);
}

/// The call `quott_mut(iox, y)`
/// stores into `iox` the truncated quotient of `iox` and `y`.
__attribute__ ((__nonnull__))
inline void $(quott_mut, A, D)(A *restrict const iox,
    A const *restrict const y) {
  for (size_t i = 0; i < D; ++i)
    iox[i] = $(quott, A)(iox[i], y[i]);
}

/// The call `remt_mut(iox, y)`
/// stores into `iox` the truncated remainder of `iox` and `y`.
/// This is analogous to the binary operator `%=`.
__attribute__ ((__nonnull__))
inline void $(remt_mut, A, D)(A *restrict const iox,
    A const *restrict const y) {
  for (size_t i = 0; i < D; ++i)
    iox[i] = $(remt, A)(iox[i], y[i]);
}

/// The call `quote_mut(iox, y)`
/// stores into `iox` the Euclidean quotient of `iox` and `y`.
/// This is analogous to the binary operator `/=`.
__attribute__ ((__nonnull__))
inline void $(quote_mut, A, D)(A *restrict const iox,
    A const *restrict const y) {
  for (size_t i = 0; i < D; ++i)
    iox[i] = $(quote, A)(iox[i], y[i]);
}

/// The call `reme_mut(iox, y)`
/// stores into `iox` the Euclidean remainder of `iox` and `y`.
/// This is analogous to the binary operator `%=`.
__attribute__ ((__nonnull__))
inline void $(reme_mut, A, D)(A *restrict const iox,
    A const *restrict const y) {
  for (size_t i = 0; i < D; ++i)
    iox[i] = $(reme, A)(iox[i], y[i]);
}
