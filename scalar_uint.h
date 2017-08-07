/// Scalar arithmetic for unsigned integer types.

/// The call `add(x, y)`
/// returns the sum of `x` and `y`.
/// This is analogous to the binary operator `+`.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A type(add, A)(A const x, A const y) {
  dynamic_assert(!(type(maxval, A)() - y < x), "Arithmetic overflow");

  return x + y;
}

/// The call `sub(x, y)`
/// returns the difference of `x` and `y`.
/// This is analogous to the binary operator `-`.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A type(sub, A)(A const x, A const y) {
  dynamic_assert(!(type(minval, A)() + y > x), "Arithmetic overflow");

  return x - y;
}

/// The call `mul(x, y)`
/// returns the product of `x` and `y`.
/// This is analogous to the binary operator `*`.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A type(mul, A)(A const x, A const y) {
  dynamic_assert(!(type(maxval, A)() / y < x), "Arithmetic overflow");

  return x * y;
}

/// The call `quot(x, y)`
/// returns the quotient of `x` and `y`.
/// This is analogous to the binary operator `/`.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A type(quot, A)(A const x, A const y) {
  dynamic_assert(!(y == type(zero, A)()), "Division by zero");

  return x / y;
}

/// The call `rem(x, y)`
/// returns the remainder of `x` and `y`.
/// This is analogous to the binary operator `%`.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A type(rem, A)(A const x, A const y) {
  dynamic_assert(!(y == type(zero, A)()), "Division by zero");

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

/// The call `quot_mut(iox, y)`
/// stores into `iox` the quotient of `iox` and `y`.
/// This is analogous to the binary operator `/=`.
__attribute__ ((__nonnull__))
inline void type(quot_mut, A)(A *const iox, A const y) {
  *iox = type(quot, A)(*iox, y);
}

/// The call `rem_mut(iox, y)`
/// stores into `iox` the remainder of `iox` and `y`.
/// This is analogous to the binary operator `%=`.
__attribute__ ((__nonnull__))
inline void type(rem_mut, A)(A *const iox, A const y) {
  *iox = type(rem, A)(*iox, y);
}
