/// Scalar arithmetic for floating-point types.

/// The call `add(x, y)`
/// returns the sum of `x` and `y`.
/// This is analogous to the binary operator `+`.
__attribute__ ((__const__, __pure__))
inline A type(add, A)(A const x, A const y) {
  return x + y;
}

/// The call `neg(x)`
/// returns the negation of `x`.
/// This is analogous to the unary operator `-`.
__attribute__ ((__const__, __pure__))
inline A type(neg, A)(A const x) {
  return -x;
}

/// The call `sub(x, y)`
/// returns the difference of `x` and `y`.
/// This is analogous to the binary operator `-`.
__attribute__ ((__const__, __pure__))
inline A type(sub, A)(A const x, A const y) {
  return x - y;
}

/// The call `mul(x, y)`
/// returns the product of `x` and `y`.
/// This is analogous to the binary operator `*`.
__attribute__ ((__const__, __pure__))
inline A type(mul, A)(A const x, A const y) {
  return x * y;
}

/// The call `recip(x)`
/// returns the reciprocal of `x`.
__attribute__ ((__const__, __pure__))
inline A type(recip, A)(A const x) {
  return type(one, A)() / x;
}

/// The call `div(x, y)`
/// returns the division of `x` and `y`.
/// This is analogous to the binary operator `/`.
__attribute__ ((__const__, __pure__))
inline A type(div, A)(A const x, A const y) {
  return x / y;
}

/// The call `quote(x, y)`
/// returns the Euclidean quotient of `x` and `y`.
__attribute__ ((__const__, __pure__))
inline A type(quote, A)(A const x, A const y) {
  A const q = type(quott, A)(x, y);
  A const r = type(remt, A)(x, y);
  A const s = r >= type(zero, A)() ?
    type(zero, A)() : y < type(zero, A)() ?
    -type(one, A)() : type(one, A)();

  return q - s;
}

/// The call `reme(x, y)`
/// returns the Euclidean remainder of `x` and `y`.
__attribute__ ((__const__, __pure__))
inline A type(reme, A)(A const x, A const y) {
  A const r = type(remt, A)(x, y);
  A const s = r >= type(zero, A)() ?
    type(zero, A)() : y < type(zero, A)() ?
    -type(one, A)() : type(one, A)();

  return r + s * y;
}

/// The call `add_mut(iox, y)`
/// stores into `iox` the sum of `iox` and `y`.
/// This is analogous to the binary operator `+=`.
__attribute__ ((__nonnull__))
inline void type(add_mut, A)(A *const iox, A const y) {
  *iox = type(add, A)(*iox, y);
}

/// The call `neg_mut(iox)`
/// stores into `iox` the negation of `iox`.
__attribute__ ((__nonnull__))
inline void type(neg_mut, A)(A *const iox) {
  *iox = type(neg, A)(*iox);
}

/// The call `sub_mut(iox, y)`
/// stores into `iox` the difference of `iox` and `y`.
/// This is analogous to the binary operator `-=`.
__attribute__ ((__nonnull__))
inline void type(sub_mut, A)(A *const iox, A const y) {
  *iox = type(sub, A)(*iox, y);
}

/// The call `mul_mut(iox, y)`
/// stores into `iox` the product of `iox` and `y`.
/// This is analogous to the binary operator `*=`.
__attribute__ ((__nonnull__))
inline void type(mul_mut, A)(A *const iox, A const y) {
  *iox = type(mul, A)(*iox, y);
}

/// The call `recip_mut(iox)`
/// stores into `iox` the reciprocal of `iox`.
__attribute__ ((__nonnull__))
inline void type(recip_mut, A)(A *const iox) {
  *iox = type(recip, A)(*iox);
}

/// The call `div_mut(iox, y)`
/// stores into `iox` the division of `iox` and `y`.
/// This is analogous to the binary operator `/=`.
__attribute__ ((__nonnull__))
inline void type(div_mut, A)(A *const iox, A const y) {
  *iox = type(div, A)(*iox, y);
}

/// The call `quott_mut(iox, y)`
/// stores into `iox` the truncated quotient of `iox` and `y`.
__attribute__ ((__nonnull__))
inline void type(quott_mut, A)(A *const iox, A const y) {
  *iox = type(quott, A)(*iox, y);
}

/// The call `remt_mut(iox, y)`
/// stores into `iox` the truncated remainder of `iox` and `y`.
__attribute__ ((__nonnull__))
inline void type(remt_mut, A)(A *const iox, A const y) {
  *iox = type(remt, A)(*iox, y);
}

/// The call `quote_mut(iox, y)`
/// stores into `iox` the Euclidean quotient of `iox` and `y`.
__attribute__ ((__nonnull__))
inline void type(quote_mut, A)(A *const iox, A const y) {
  *iox = type(quote, A)(*iox, y);
}

/// The call `reme_mut(iox, y)`
/// stores into `iox` the Euclidean remainder of `iox` and `y`.
__attribute__ ((__nonnull__))
inline void type(reme_mut, A)(A *const iox, A const y) {
  *iox = type(reme, A)(*iox, y);
}
