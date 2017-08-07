/// Scalar arithmetic for integer types.

/// The call `succ(x)`
/// returns the successor of `x`.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A $(succ, A)(A const x) {
  dynamic_assert(!(x == $(maxval, A)()), "Arithmetic overflow");

  return x + $(one, A)();
}

/// The call `pred(x)`
/// returns the predecessor of `x`.
#ifndef DEBUG
__attribute__ ((__const__, __pure__))
#endif
inline A $(pred, A)(A const x) {
  dynamic_assert(!(x == $(minval, A)()), "Arithmetic overflow");

  return x - $(one, A)();
}

/// The call `succ_mut(iox)`
/// stores into `iox` the successor of `iox`.
/// This is analogous to the unary operator `++`.
__attribute__ ((__nonnull__))
inline void $(succ_mut, A)(A *const iox) {
  *iox = $(succ, A)(*iox);
}

/// The call `pred_mut(iox)`
/// stores into `iox` the predecessor of `iox`.
/// This is analogous to the unary operator `--`.
__attribute__ ((__nonnull__))
inline void $(pred_mut, A)(A *const iox) {
  *iox = $(pred, A)(*iox);
}
