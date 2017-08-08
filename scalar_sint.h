/// Scalar arithmetic for signed integer types.

/// The call `add_ovf(x, y)`
/// checks whether the sum of `x` and `y` would overflow.
__attribute__ ((__const__, __pure__))
inline bool $(add_ovf, A)(A const x, A const y) {
  return (y > $(zero, A)() && $(maxval, A)() - y < x) ||
    (y < $(zero, A)() && $(minval, A)() - y > x);
}

/// The call `neg_ovf(x)`
/// checks whether the negation of `x` would overflow.
__attribute__ ((__const__, __pure__))
inline bool $(neg_ovf, A)(A const x) {
  return (x > $(zero, A)() && $(minval, A)() + x > $(zero, A)()) ||
    (x < $(zero, A)() && $(maxval, A)() + x < $(zero, A)());
}

/// The call `sub_ovf(x, y)`
/// checks whether the difference of `x` and `y` would overflow.
__attribute__ ((__const__, __pure__))
inline bool $(sub_ovf, A)(A const x, A const y) {
  return (y > $(zero, A)() && $(minval, A)() + y > x) ||
    (y < $(zero, A)() && $(maxval, A)() + y < x);
}

/// The call `mul_ovf(x, y)`
/// checks whether the product of `x` and `y` would overflow.
__attribute__ ((__const__, __pure__))
inline bool $(mul_ovf, A)(A const x, A const y) {
  return x != $(zero, A)() && (x > $(zero, A)() ?
      (y > $(zero, A)() ? $(maxval, A)() / x < y :
        $(minval, A)() / x > y) :
      (y > $(zero, A)() ? $(minval, A)() / y > x :
        ($(minval, A)() + $(maxval, A)() <= $(zero, A)() ?
          $(maxval, A)() / x > y :
          $(maxval, A)() / -x < -y)));
}

/// The call `quott_ovf(x, y)`
/// checks whether the truncated quotient of `x` and `y` would overflow.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline bool $(quott_ovf, A)(A const x, A const y) {
  dynamic_assert(
      ($(minval, A)() + $(maxval, A)() <= $(zero, A)() &&
        $(minval, A)() / ($(one, A)() + $(one, A)()) +
        $(maxval, A)() >= $(zero, A)()) ||
      ($(minval, A)() + $(maxval, A)() >= $(zero, A)() &&
        $(minval, A)() +
        $(maxval, A)() / ($(one, A)() + $(one, A)()) <= $(zero, A)()),
      "Unsupported limits");

  return $(minval, A)() + $(maxval, A)() <= $(zero, A)() ?
    x < -$(maxval, A)() && y < $(zero, A)() && y == -$(one, A)() :
    x > -$(minval, A)() && y < $(zero, A)() && y == -$(one, A)();
}

/// The call `remt_ovf(x, y)`
/// checks whether the truncated remainder of `x` and `y` would overflow.
__attribute__ ((__const__, __pure__))
inline bool $(remt_ovf, A)(__attribute__ ((__unused__)) A const x,
    __attribute__ ((__unused__)) A const y) {
  return false;
}

/// The call `add(x, y)`
/// returns the sum of `x` and `y`.
/// This is analogous to the binary operator `+`.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A $(add, A)(A const x, A const y) {
  dynamic_assert(!$(add_ovf, A)(x, y), "Arithmetic overflow");

  return x + y;
}

/// The call `neg(x)`
/// returns the negation of `x`.
/// This is analogous to the unary operator `-`.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A $(neg, A)(A const x) {
  dynamic_assert(!$(neg_ovf, A)(x), "Arithmetic overflow");

  return -x;
}

/// The call `sub(x, y)`
/// returns the difference of `x` and `y`.
/// This is analogous to the binary operator `-`.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A $(sub, A)(A const x, A const y) {
  dynamic_assert(!$(sub_ovf, A)(x, y), "Arithmetic overflow");

  return x - y;
}

/// The call `mul(x, y)`
/// returns the product of `x` and `y`.
/// This is analogous to the binary operator `*`.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A $(mul, A)(A const x, A const y) {
  dynamic_assert(!$(mul_ovf, A)(x, y), "Arithmetic overflow");

  return x * y;
}

/// The call `quott(x, y)`
/// returns the truncated quotient of `x` and `y`.
/// This is analogous to the binary operator `/`.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A $(quott, A)(A const x, A const y) {
  dynamic_assert(y != $(zero, A)(), "Division by zero");
  dynamic_assert(!$(quott_ovf, A)(x, y), "Arithmetic overflow");

  return x / y;
}

/// The call `remt(x, y)`
/// returns the truncated remainder of `x` and `y`.
/// This is analogous to the binary operator `%`.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A $(remt, A)(A const x, A const y) {
  dynamic_assert(y != $(zero, A)(), "Division by zero");
  dynamic_assert(!$(remt_ovf, A)(x, y), "Arithmetic overflow");

  return x % y;
}

/// The call `quotremt(oz, x, y)`
/// stores into `oz` the truncated quotient and remainder of `x` and `y`.
__attribute__ ((__nonnull__))
inline void $(quotremt, A)(A *const oz, A const x, A const y) {
  oz[0] = $(quott, A)(x, y);
  oz[1] = $(remt, A)(x, y);
}

/// The call `quote(x, y)`
/// returns the Euclidean quotient of `x` and `y`.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A $(quote, A)(A const x, A const y) {
  A const q = $(quott, A)(x, y);
  A const r = $(remt, A)(x, y);
  A const s = r >= $(zero, A)() ?
    $(zero, A)() : y < $(zero, A)() ?
    -$(one, A)() : $(one, A)();

  return $(sub, A)(q, s);
}

/// The call `reme(x, y)`
/// returns the Euclidean remainder of `x` and `y`.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A $(reme, A)(A const x, A const y) {
  A const r = $(remt, A)(x, y);
  A const s = r >= $(zero, A)() ?
    $(zero, A)() : y < $(zero, A)() ?
    -$(one, A)() : $(one, A)();

  return $(add, A)(r, $(mul, A)(s, y));
}

/// The call `quotreme(oz, x, y)`
/// stores into `oz` the Euclidean quotient and remainder of `x` and `y`.
__attribute__ ((__nonnull__))
inline void $(quotreme, A)(A *const oz, A const x, A const y) {
  A const q = $(quott, A)(x, y);
  A const r = $(remt, A)(x, y);
  A const s = r >= $(zero, A)() ?
    $(zero, A)() : y < $(zero, A)() ?
    -$(one, A)() : $(one, A)();

  oz[0] = $(sub, A)(q, s);
  oz[1] = $(add, A)(r, $(mul, A)(s, y));
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

/// The call `quott_mut(iox, y)`
/// stores into `iox` the truncated quotient of `iox` and `y`.
/// This is analogous to the binary operator `/=`.
__attribute__ ((__nonnull__))
inline void $(quott_mut, A)(A *const iox, A const y) {
  *iox = $(quott, A)(*iox, y);
}

/// The call `remt_mut(iox, y)`
/// stores into `iox` the truncated remainder of `iox` and `y`.
/// This is analogous to the binary operator `%=`.
__attribute__ ((__nonnull__))
inline void $(remt_mut, A)(A *const iox, A const y) {
  *iox = $(remt, A)(*iox, y);
}

/// The call `quotremt_mut(iox)`
/// stores into `iox` the truncated quotient and remainder of `iox`.
__attribute__ ((__nonnull__))
inline void $(quotremt_mut, A)(A *const iox) {
  $(quotremt, A)(iox, iox[0], iox[1]);
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

/// The call `quotreme_mut(iox)`
/// stores into `iox` the truncated quotient and remainder of `iox`.
__attribute__ ((__nonnull__))
inline void $(quotreme_mut, A)(A *const iox) {
  $(quotreme, A)(iox, iox[0], iox[1]);
}
