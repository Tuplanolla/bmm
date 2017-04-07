/// Preprocessor directives.
#ifndef BMM_CPP_H
#define BMM_CPP_H

/// The preprocessor directive `MIN(x, y)`
/// expands to the lesser of `x` and `y`.
#define MIN(x, y) ((x) < (y) ? (x) : (y))

/// The preprocessor directive `MAX(x, y)`
/// expands to the greater of `x` and `y`.
#define MAX(x, y) ((x) > (y) ? (x) : (y))

#endif
