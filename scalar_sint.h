/// Scalar arithmetic for signed integer types.

/// The call `add(x, y)`
/// returns the sum of `x` and `y`.
/// This is analogous to the binary operator `+`.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A type(add, A)(A const x, A const y) {
  dynamic_assert(
      !((y > type(zero, A)() && type(maxval, A)() - y < x) ||
        (y < type(zero, A)() && type(minval, A)() - y > x)),
      "Arithmetic overflow");

  return x + y;
}

/// The call `neg(x)`
/// returns the negation of `x`.
/// This is analogous to the unary operator `-`.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A type(neg, A)(A const x) {
  dynamic_assert(
      !((x > type(zero, A)() && type(minval, A)() + x > type(zero, A)()) ||
        (x < type(zero, A)() && type(maxval, A)() + x < type(zero, A)())),
      "Arithmetic overflow");

  return -x;
}

/// The call `sub(x, y)`
/// returns the difference of `x` and `y`.
/// This is analogous to the binary operator `-`.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A type(sub, A)(A const x, A const y) {
  dynamic_assert(
      !((y > type(zero, A)() && type(minval, A)() + y > x) ||
        (y < type(zero, A)() && type(maxval, A)() + y < x)),
      "Arithmetic overflow");

  return x - y;
}

// It is very difficult to reliably detect overflows here,
// because the standard does not specify
// whether negative limits have larger magnitudes
// than the corresponding positive limits or vice versa.
// This cannot be easily computationally deduced either,
// since comparing the magnitudes of the limits without an overflow
// would first require negating the one with a smaller magnitude.
// For this reason we assume that for multiplication and division
// `type(minval, A)() <= -type(maxval, A)()` and for division
// `type(minval, A)() / type(two, A)() >= -type(maxval, A)())`.

/// The call `mul(x, y)`
/// returns the product of `x` and `y`.
/// This is analogous to the binary operator `*`.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A type(mul, A)(A const x, A const y) {
  dynamic_assert(!(y == type(zero, A)()), "Division by zero");
  dynamic_assert(
      !(x > type(zero, A)() ? (y > type(zero, A)() ?
          type(maxval, A)() / y < x : type(minval, A)() / x > y) :
        x < type(zero, A)() ? (y > type(zero, A)() ?
          type(minval, A)() / y > x : type(maxval, A)() / x > y) : false),
      "Arithmetic overflow");

  return x - y;
}

/// The call `quot(x, y)`
/// returns the quotient of `x` and `y`.
/// This is analogous to the binary operator `/`.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A type(quot, A)(A const x, A const y) {
  dynamic_assert(!(y == type(zero, A)()), "Division by zero");
  dynamic_assert(!(x < -type(maxval, A)() &&
        y < type(zero, A)() && y == -type(one, A)()),
      "Arithmetic overflow");

  A const q = x / y;
  A const r = x % y;
  A const s = r >= type(zero, A)() ?
    type(zero, A)() : y < type(zero, A)() ?
    -type(one, A)() : type(one, A)();

  return q - s;
}

/// The call `rem(x, y)`
/// returns the remainder of `x` and `y`.
/// This is analogous to the binary operator `%`.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A type(rem, A)(A const x, A const y) {
  dynamic_assert(!(y == type(zero, A)()), "Division by zero");
  dynamic_assert(!(x < -type(maxval, A)() &&
        y < type(zero, A)() && y == -type(one, A)()),
      "Arithmetic overflow");

  A const r = x % y;
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
