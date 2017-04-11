#ifndef BMM_CPP_H
/// Preprocessor directives.
#define BMM_CPP_H

/// The preprocessor directive `BMM_MIN(x, y)`
/// expands to the lesser of `x` and `y`.
#define BMM_MIN(x, y) ((x) < (y) ? (x) : (y))

/// The preprocessor directive `BMM_MAX(x, y)`
/// expands to the greater of `x` and `y`.
#define BMM_MAX(x, y) ((x) > (y) ? (x) : (y))

/// The preprocessor directive `BMM_TESTBIT(x, p)`
/// checks whether `x` has the bit at exponent `p` set.
#define BMM_TESTBIT(x, p) ((((x) >> (p)) & 1) != 0)

/// The preprocessor directive `BMM_SETBIT(x, p)`
/// expands to `x` with the bit at exponent `p` set.
#define BMM_SETBIT(x, p) ((x) | (1 << (p)))

/// The preprocessor directive `BMM_CLEARBIT(x, p)`
/// expands to `x` with the bit at exponent `p` unset.
#define BMM_CLEARBIT(x, p) ((x) & ~(1 << (p)))

#endif
