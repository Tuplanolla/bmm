/// Scalar arithmetic for unsigned integer types.

/// The call `add_ovf(x, y)`
/// checks whether the sum of `x` and `y` would overflow.
__attribute__ ((__const__, __pure__))
inline bool type(add_ovf, A)(A const x, A const y) {
  return type(maxval, A)() - y < x;
}

/// The call `sub_ovf(x, y)`
/// checks whether the difference of `x` and `y` would overflow.
__attribute__ ((__const__, __pure__))
inline bool type(sub_ovf, A)(A const x, A const y) {
  return y > x;
}

/// The call `mul_ovf(x, y)`
/// checks whether the product of `x` and `y` would overflow.
__attribute__ ((__const__, __pure__))
inline bool type(mul_ovf, A)(A const x, A const y) {
  return x != type(zero, A)() && type(maxval, A)() / x < y;
}

/// The call `quott_ovf(x, y)`
/// checks whether the truncated quotient of `x` and `y` would overflow.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline bool type(quott_ovf, A)(__attribute__ ((__unused__)) A const x,
    __attribute__ ((__unused__)) A const y) {
  return false;
}

/// The call `remt_ovf(x, y)`
/// checks whether the truncated remainder of `x` and `y` would overflow.
__attribute__ ((__const__, __pure__))
inline bool type(remt_ovf, A)(__attribute__ ((__unused__)) A const x,
    __attribute__ ((__unused__)) A const y) {
  return false;
}

/// The call `quote_ovf(x, y)`
/// checks whether the Euclidean quotient of `x` and `y` would overflow.
__attribute__ ((__const__, __pure__))
inline bool type(quote_ovf, A)(__attribute__ ((__unused__)) A const x,
    __attribute__ ((__unused__)) A const y) {
  return false;
}

/// The call `reme_ovf(x, y)`
/// checks whether the Euclidean remainder of `x` and `y` would overflow.
__attribute__ ((__const__, __pure__))
inline bool type(reme_ovf, A)(__attribute__ ((__unused__)) A const x,
    __attribute__ ((__unused__)) A const y) {
  return false;
}

/// The call `add(x, y)`
/// returns the sum of `x` and `y`.
/// This is analogous to the binary operator `+`.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A type(add, A)(A const x, A const y) {
  dynamic_assert(!type(add_ovf, A)(x, y), "Arithmetic overflow");

  return x + y;
}

/// The call `sub(x, y)`
/// returns the difference of `x` and `y`.
/// This is analogous to the binary operator `-`.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A type(sub, A)(A const x, A const y) {
  dynamic_assert(!type(sub_ovf, A)(x, y), "Arithmetic overflow");

  return x - y;
}

/// The call `mul(x, y)`
/// returns the product of `x` and `y`.
/// This is analogous to the binary operator `*`.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A type(mul, A)(A const x, A const y) {
  dynamic_assert(!type(mul_ovf, A)(x, y), "Arithmetic overflow");

  return x * y;
}

/// The call `quott(x, y)`
/// returns the truncated quotient of `x` and `y`.
/// This is analogous to the binary operator `/`.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A type(quott, A)(A const x, A const y) {
  dynamic_assert(y != type(zero, A)(), "Division by zero");
  dynamic_assert(!type(quott_ovf, A)(x, y), "Arithmetic overflow");

  return x / y;
}

/// The call `remt(x, y)`
/// returns the truncated remainder of `x` and `y`.
/// This is analogous to the binary operator `%`.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A type(remt, A)(A const x, A const y) {
  dynamic_assert(y != type(zero, A)(), "Division by zero");
  dynamic_assert(!type(remt_ovf, A)(x, y), "Arithmetic overflow");

  return x % y;
}

/// The call `quote(x, y)`
/// returns the Euclidean quotient of `x` and `y`.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A type(quote, A)(A const x, A const y) {
  dynamic_assert(y != type(zero, A)(), "Division by zero");
  dynamic_assert(!type(quote_ovf, A)(x, y), "Arithmetic overflow");

  return x / y;
}

/// The call `reme(x, y)`
/// returns the Euclidean remainder of `x` and `y`.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A type(reme, A)(A const x, A const y) {
  dynamic_assert(y != type(zero, A)(), "Division by zero");
  dynamic_assert(!type(reme_ovf, A)(x, y), "Arithmetic overflow");

  return x % y;
}

/// The call `add_mut(iox, y)`
/// stores into `iox` the sum of `iox` and `y`.
/// This is analogous to the binary operator `+=`.
__attribute__ ((__nonnull__))
inline void type(add_mut, A)(A *const iox, A const y) {
  *iox = type(add, A)(*iox, y);
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

/// The call `quott_mut(iox, y)`
/// stores into `iox` the truncated quotient of `iox` and `y`.
/// This is analogous to the binary operator `/=`.
__attribute__ ((__nonnull__))
inline void type(quott_mut, A)(A *const iox, A const y) {
  *iox = type(quott, A)(*iox, y);
}

/// The call `remt_mut(iox, y)`
/// stores into `iox` the truncated remainder of `iox` and `y`.
/// This is analogous to the binary operator `%=`.
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
