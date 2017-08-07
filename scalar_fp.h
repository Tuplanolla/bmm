/// Scalar arithmetic for floating-point types.

/// The call `add(x, y)`
/// returns the sum of `x` and `y`.
/// This is analogous to the binary operator `+`.
__attribute__ ((__const__, __pure__))
inline A $(add, A)(A const x, A const y) {
  return x + y;
}

/// The call `neg(x)`
/// returns the negation of `x`.
/// This is analogous to the unary operator `-`.
__attribute__ ((__const__, __pure__))
inline A $(neg, A)(A const x) {
  return -x;
}

/// The call `sub(x, y)`
/// returns the difference of `x` and `y`.
/// This is analogous to the binary operator `-`.
__attribute__ ((__const__, __pure__))
inline A $(sub, A)(A const x, A const y) {
  return x - y;
}

/// The call `mul(x, y)`
/// returns the product of `x` and `y`.
/// This is analogous to the binary operator `*`.
__attribute__ ((__const__, __pure__))
inline A $(mul, A)(A const x, A const y) {
  return x * y;
}

/// The call `recip(x)`
/// returns the reciprocal of `x`.
__attribute__ ((__const__, __pure__))
inline A $(recip, A)(A const x) {
  return $(one, A)() / x;
}

/// The call `div(x, y)`
/// returns the division of `x` and `y`.
/// This is analogous to the binary operator `/`.
__attribute__ ((__const__, __pure__))
inline A $(div, A)(A const x, A const y) {
  return x / y;
}

/// The call `quote(x, y)`
/// returns the Euclidean quotient of `x` and `y`.
__attribute__ ((__const__, __pure__))
inline A $(quote, A)(A const x, A const y) {
  A const q = $(quott, A)(x, y);
  A const r = $(remt, A)(x, y);
  A const s = r >= $(zero, A)() ?
    $(zero, A)() : y < $(zero, A)() ?
    -$(one, A)() : $(one, A)();

  return q - s;
}

/// The call `reme(x, y)`
/// returns the Euclidean remainder of `x` and `y`.
__attribute__ ((__const__, __pure__))
inline A $(reme, A)(A const x, A const y) {
  A const r = $(remt, A)(x, y);
  A const s = r >= $(zero, A)() ?
    $(zero, A)() : y < $(zero, A)() ?
    -$(one, A)() : $(one, A)();

  return r + s * y;
}

/// The call `add_mut(iox, y)`
/// stores into `iox` the sum of `iox` and `y`.
/// This is analogous to the binary operator `+=`.
__attribute__ ((__nonnull__))
inline void $(add_mut, A)(A *const iox, A const y) {
  *iox = $(add, A)(*iox, y);
}

/// The call `neg_mut(iox)`
/// stores into `iox` the negation of `iox`.
__attribute__ ((__nonnull__))
inline void $(neg_mut, A)(A *const iox) {
  *iox = $(neg, A)(*iox);
}

/// The call `sub_mut(iox, y)`
/// stores into `iox` the difference of `iox` and `y`.
/// This is analogous to the binary operator `-=`.
__attribute__ ((__nonnull__))
inline void $(sub_mut, A)(A *const iox, A const y) {
  *iox = $(sub, A)(*iox, y);
}

/// The call `mul_mut(iox, y)`
/// stores into `iox` the product of `iox` and `y`.
/// This is analogous to the binary operator `*=`.
__attribute__ ((__nonnull__))
inline void $(mul_mut, A)(A *const iox, A const y) {
  *iox = $(mul, A)(*iox, y);
}

/// The call `recip_mut(iox)`
/// stores into `iox` the reciprocal of `iox`.
__attribute__ ((__nonnull__))
inline void $(recip_mut, A)(A *const iox) {
  *iox = $(recip, A)(*iox);
}

/// The call `div_mut(iox, y)`
/// stores into `iox` the division of `iox` and `y`.
/// This is analogous to the binary operator `/=`.
__attribute__ ((__nonnull__))
inline void $(div_mut, A)(A *const iox, A const y) {
  *iox = $(div, A)(*iox, y);
}

/// The call `quott_mut(iox, y)`
/// stores into `iox` the truncated quotient of `iox` and `y`.
__attribute__ ((__nonnull__))
inline void $(quott_mut, A)(A *const iox, A const y) {
  *iox = $(quott, A)(*iox, y);
}

/// The call `remt_mut(iox, y)`
/// stores into `iox` the truncated remainder of `iox` and `y`.
__attribute__ ((__nonnull__))
inline void $(remt_mut, A)(A *const iox, A const y) {
  *iox = $(remt, A)(*iox, y);
}

/// The call `quote_mut(iox, y)`
/// stores into `iox` the Euclidean quotient of `iox` and `y`.
__attribute__ ((__nonnull__))
inline void $(quote_mut, A)(A *const iox, A const y) {
  *iox = $(quote, A)(*iox, y);
}

/// The call `reme_mut(iox, y)`
/// stores into `iox` the Euclidean remainder of `iox` and `y`.
__attribute__ ((__nonnull__))
inline void $(reme_mut, A)(A *const iox, A const y) {
  *iox = $(reme, A)(*iox, y);
}
