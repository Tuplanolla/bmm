/// Scalar arithmetic for lifted types.

/// The call `zero()`
/// returns zero.
/// This is analogous to the constant `0`.
__attribute__ ((__const__, __pure__))
inline signed_char type(zero, signed_char)(void) {
  return 0;
}

__attribute__ ((__const__, __pure__))
inline unsigned_char type(zero, unsigned_char)(void) {
  return 0;
}

__attribute__ ((__const__, __pure__))
inline int type(zero, int)(void) {
  return 0;
}

__attribute__ ((__const__, __pure__))
inline double type(zero, double)(void) {
  return 0.0;
}

__attribute__ ((__const__, __pure__))
inline size_t type(zero, size_t)(void) {
  return 0;
}

/// The call `one()`
/// returns one.
/// This is analogous to the constant `1`.
__attribute__ ((__const__, __pure__))
inline signed_char type(one, signed_char)(void) {
  return 1;
}

__attribute__ ((__const__, __pure__))
inline unsigned_char type(one, unsigned_char)(void) {
  return 1;
}

__attribute__ ((__const__, __pure__))
inline int type(one, int)(void) {
  return 1;
}

__attribute__ ((__const__, __pure__))
inline double type(one, double)(void) {
  return 1.0;
}

__attribute__ ((__const__, __pure__))
inline size_t type(one, size_t)(void) {
  return 1;
}

/// The call `minval()`
/// returns the minimal representable value.
/// This is analogous to the constant `INT_MIN` for `int`.
__attribute__ ((__const__, __pure__))
inline signed_char type(minval, signed_char)(void) {
  return SCHAR_MIN;
}

__attribute__ ((__const__, __pure__))
inline unsigned_char type(minval, unsigned_char)(void) {
  return 0;
}

__attribute__ ((__const__, __pure__))
inline int type(minval, int)(void) {
  return INT_MIN;
}

__attribute__ ((__const__, __pure__))
inline double type(minval, double)(void) {
  return -DBL_MAX;
}

__attribute__ ((__const__, __pure__))
inline size_t type(minval, size_t)(void) {
  return 0;
}

/// The call `maxval()`
/// returns the maximal representable value.
/// This is analogous to the constant `INT_MAX` for `int`.
__attribute__ ((__const__, __pure__))
inline signed_char type(maxval, signed_char)(void) {
  return SCHAR_MAX;
}

__attribute__ ((__const__, __pure__))
inline unsigned_char type(maxval, unsigned_char)(void) {
  return UCHAR_MAX;
}

__attribute__ ((__const__, __pure__))
inline int type(maxval, int)(void) {
  return INT_MAX;
}

__attribute__ ((__const__, __pure__))
inline double type(maxval, double)(void) {
  return DBL_MAX;
}

__attribute__ ((__const__, __pure__))
inline size_t type(maxval, size_t)(void) {
  return SIZE_MAX;
}

/// The call `bmm_trunc(x, y)`
/// returns the truncated quotient of `x` and `y`.
__attribute__ ((__const__, __pure__))
inline signed_char type(trunc, signed_char)(signed_char const x,
    signed_char const y) {
  return x / y;
}

__attribute__ ((__const__, __pure__))
inline unsigned_char type(trunc, unsigned_char)(unsigned_char const x,
    unsigned_char const y) {
  return x / y;
}

__attribute__ ((__const__, __pure__))
inline int type(trunc, int)(int const x, int const y) {
  return x / y;
}

__attribute__ ((__const__, __pure__))
inline double type(trunc, double)(double const x, double const y) {
  return trunc(x / y);
}

__attribute__ ((__const__, __pure__))
inline size_t type(trunc, size_t)(size_t const x, size_t const y) {
  return x / y;
}

/// The call `bmm_mod(x)`
/// returns the truncated remainder of `x` and `y`.
__attribute__ ((__const__, __pure__))
inline signed_char type(mod, signed_char)(signed_char const x,
    signed_char const y) {
  return x % y;
}

__attribute__ ((__const__, __pure__))
inline unsigned_char type(mod, unsigned_char)(unsigned_char const x,
    unsigned_char const y) {
  return x % y;
}

__attribute__ ((__const__, __pure__))
inline int type(mod, int)(int const x, int const y) {
  return x % y;
}

__attribute__ ((__const__, __pure__))
inline double type(mod, double)(double const x, double const y) {
  return fmod(x, y);
}

__attribute__ ((__const__, __pure__))
inline size_t type(mod, size_t)(size_t const x, size_t const y) {
  return x % y;
}
