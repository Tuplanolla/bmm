/// The call `zero(ox)`
/// stores into `ox` zero.
/// This is analogous to the constant `0`.
__attribute__ ((__nonnull__))
inline void type(zero, A, D)(A *const ox) {
  for (size_t i = 0; i < D; ++i)
    ox[i] = type(zero, A)();
}

/// The call `add(oz, x, y)`
/// stores into `oz` the sum of `x` and `y`.
/// This is analogous to the binary operator `+`.
__attribute__ ((__nonnull__))
inline void type(add, A, D)(A *restrict const oz,
    A const *restrict const x, A const *restrict const y) {
  for (size_t i = 0; i < D; ++i)
    oz[i] = type(add, A)(x[i], y[i]);
}

/// The call `neg(oy, x)`
/// stores into `oy` the negation of `x`.
/// This is analogous to the unary operator `-`.
__attribute__ ((__nonnull__))
inline void type(neg, A, D)(A *restrict const oy, A const *restrict const x) {
  for (size_t i = 0; i < D; ++i)
    oy[i] = type(neg, A)(x[i]);
}

/// The call `sub(oz, x, y)`
/// stores into `oz` the difference of `x` and `y`.
/// This is analogous to the binary operator `-`.
__attribute__ ((__nonnull__))
inline void type(sub, A, D)(A *restrict const oz,
    A const *restrict const x, A const *restrict const y) {
  for (size_t i = 0; i < D; ++i)
    oz[i] = type(sub, A)(x[i], y[i]);
}

/// The call `one(ox)`
/// stores into `ox` one.
/// This is analogous to the constant `1`.
__attribute__ ((__nonnull__))
inline void type(one, A, D)(A *const ox) {
  for (size_t i = 0; i < D; ++i)
    ox[i] = type(one, A)();
}

/// The call `mul(oz, x, y)`
/// stores into `oz` the product of `x` and `y`.
/// This is analogous to the binary operator `*`.
__attribute__ ((__nonnull__))
inline void type(mul, A, D)(A *restrict const oz,
    A const *restrict const x, A const *restrict const y) {
  for (size_t i = 0; i < D; ++i)
    oz[i] = type(mul, A)(x[i], y[i]);
}

/// The call `recip(oy, x)`
/// stores into `oy` the reciprocal of `x`.
__attribute__ ((__nonnull__))
inline void type(recip, A, D)(A *restrict const oy,
    A const *restrict const x) {
  for (size_t i = 0; i < D; ++i)
    oy[i] = type(recip, A)(x[i]);
}

/// The call `div(oz, x, y)`
/// stores into `oz` the quotient of `x` and `y`.
/// This is analogous to the binary operator `/`.
__attribute__ ((__nonnull__))
inline void type(div, A, D)(A *restrict const oz,
    A const *restrict const x, A const *restrict const y) {
  for (size_t i = 0; i < D; ++i)
    oz[i] = type(div, A)(x[i], y[i]);
}

/// The call `quot(oz, x, y)`
/// stores into `oz` the quotient of `x` and `y`.
/// This is analogous to the binary operator `/`.
__attribute__ ((__nonnull__))
inline void type(quot, A, D)(A *restrict const oz,
    A const *restrict const x, A const *restrict const y) {
  for (size_t i = 0; i < D; ++i)
    oz[i] = type(quot, A)(x[i], y[i]);
}

/// The call `rem(oz, x, y)`
/// stores into `oz` the remainder of `x` and `y`.
/// This is analogous to the binary operator `%`.
__attribute__ ((__nonnull__))
inline void type(rem, A, D)(A *restrict const oz,
    A const *restrict const x, A const *restrict const y) {
  for (size_t i = 0; i < D; ++i)
    oz[i] = type(rem, A)(x[i], y[i]);
}

/// The call `add_mut(iox, y)`
/// stores into `iox` the sum of `iox` and `y`.
/// This is analogous to the binary operator `+=`.
__attribute__ ((__nonnull__))
inline void type(add_mut, A, D)(A *restrict const iox,
    A const *restrict const y) {
  for (size_t i = 0; i < D; ++i)
    iox[i] = type(add, A)(iox[i], y[i]);
}

/// The call `neg_mut(iox)`
/// stores into `iox` the negation of `iox`.
__attribute__ ((__nonnull__))
inline void type(neg_mut, A, D)(A *restrict const iox) {
  for (size_t i = 0; i < D; ++i)
    iox[i] = type(neg, A)(iox[i]);
}

/// The call `sub_mut(iox, y)`
/// stores into `iox` the difference of `iox` and `y`.
/// This is analogous to the binary operator `-=`.
__attribute__ ((__nonnull__))
inline void type(sub_mut, A, D)(A *restrict const iox,
    A const *restrict const y) {
  for (size_t i = 0; i < D; ++i)
    iox[i] = type(sub, A)(iox[i], y[i]);
}

/// The call `mul_mut(iox, y)`
/// stores into `iox` the product of `iox` and `y`.
/// This is analogous to the binary operator `*=`.
__attribute__ ((__nonnull__))
inline void type(mul_mut, A, D)(A *restrict const iox,
    A const *restrict const y) {
  for (size_t i = 0; i < D; ++i)
    iox[i] = type(mul, A)(iox[i], y[i]);
}

/// The call `recip_mut(iox)`
/// stores into `iox` the reciprocal of `iox`.
__attribute__ ((__nonnull__))
inline void type(recip_mut, A, D)(A *restrict const iox) {
  for (size_t i = 0; i < D; ++i)
    iox[i] = type(recip, A)(iox[i]);
}

/// The call `div_mut(iox, y)`
/// stores into `iox` the division of `iox` and `y`.
/// This is analogous to the binary operator `/`.
__attribute__ ((__nonnull__))
inline void type(div_mut, A, D)(A *restrict const iox,
    A const *restrict const y) {
  for (size_t i = 0; i < D; ++i)
    iox[i] = type(div, A)(iox[i], y[i]);
}

/// The call `quot_mut(iox, y)`
/// stores into `iox` the quotient of `iox` and `y`.
/// This is analogous to the binary operator `/=`.
__attribute__ ((__nonnull__))
inline void type(quot_mut, A, D)(A *restrict const iox,
    A const *restrict const y) {
  for (size_t i = 0; i < D; ++i)
    iox[i] = type(quot, A)(iox[i], y[i]);
}

/// The call `rem_mut(iox, y)`
/// stores into `iox` the remainder of `iox` and `y`.
/// This is analogous to the binary operator `%=`.
__attribute__ ((__nonnull__))
inline void type(rem_mut, A, D)(A *restrict const iox,
    A const *restrict const y) {
  for (size_t i = 0; i < D; ++i)
    iox[i] = type(rem, A)(iox[i], y[i]);
}
