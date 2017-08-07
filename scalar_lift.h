/// Scalar arithmetic for lifted types.

/// The call `zero()`
/// returns zero.
/// This is analogous to the constant `0`.
__attribute__ ((__const__, __pure__))
inline signed_char $(zero, signed_char)(void) {
  return 0;
}

__attribute__ ((__const__, __pure__))
inline unsigned_char $(zero, unsigned_char)(void) {
  return 0;
}

__attribute__ ((__const__, __pure__))
inline int $(zero, int)(void) {
  return 0;
}

__attribute__ ((__const__, __pure__))
inline double $(zero, double)(void) {
  return 0.0;
}

__attribute__ ((__const__, __pure__))
inline size_t $(zero, size_t)(void) {
  return 0;
}

/// The call `one()`
/// returns one.
/// This is analogous to the constant `1`.
__attribute__ ((__const__, __pure__))
inline signed_char $(one, signed_char)(void) {
  return 1;
}

__attribute__ ((__const__, __pure__))
inline unsigned_char $(one, unsigned_char)(void) {
  return 1;
}

__attribute__ ((__const__, __pure__))
inline int $(one, int)(void) {
  return 1;
}

__attribute__ ((__const__, __pure__))
inline double $(one, double)(void) {
  return 1.0;
}

__attribute__ ((__const__, __pure__))
inline size_t $(one, size_t)(void) {
  return 1;
}

/// The call `minval()`
/// returns the minimal representable value.
/// This is analogous to the constant `INT_MIN` for `int`.
__attribute__ ((__const__, __pure__))
inline signed_char $(minval, signed_char)(void) {
  return SCHAR_MIN;
}

__attribute__ ((__const__, __pure__))
inline unsigned_char $(minval, unsigned_char)(void) {
  return 0;
}

__attribute__ ((__const__, __pure__))
inline int $(minval, int)(void) {
  return INT_MIN;
}

__attribute__ ((__const__, __pure__))
inline double $(minval, double)(void) {
  return -DBL_MAX;
}

__attribute__ ((__const__, __pure__))
inline size_t $(minval, size_t)(void) {
  return 0;
}

/// The call `maxval()`
/// returns the maximal representable value.
/// This is analogous to the constant `INT_MAX` for `int`.
__attribute__ ((__const__, __pure__))
inline signed_char $(maxval, signed_char)(void) {
  return SCHAR_MAX;
}

__attribute__ ((__const__, __pure__))
inline unsigned_char $(maxval, unsigned_char)(void) {
  return UCHAR_MAX;
}

__attribute__ ((__const__, __pure__))
inline int $(maxval, int)(void) {
  return INT_MAX;
}

__attribute__ ((__const__, __pure__))
inline double $(maxval, double)(void) {
  return DBL_MAX;
}

__attribute__ ((__const__, __pure__))
inline size_t $(maxval, size_t)(void) {
  return SIZE_MAX;
}

/// The call `quott(x, y)`
/// returns the truncated quotient of `x` and `y`.
/// This is analogous to the binary operator `/`.
__attribute__ ((__const__, __pure__))
inline double $(quott, double)(double const x, double const y) {
  return trunc(x / y);
}

/// The call `remt(x, y)`
/// returns the truncated remainder of `x` and `y`.
/// This is analogous to the binary operator `%`.
__attribute__ ((__const__, __pure__))
inline double $(remt, double)(double const x, double const y) {
  return fmod(x, y);
}
