// Language extensions.
#ifndef BMM_EXT_H
#define BMM_EXT_H

#include <assert.h>
#include <stdbool.h>

// The preprocessor directive `ever` makes it possible to write `for ever`.
// This is essential.
#ifndef ever
#define ever (;;)
#endif

// The preprocessor directives `BEGIN` and `END` allow converting
// multiple statements into a single statement
// with local scope and no redundant semicolons.
// They must always appear together.
#define BEGIN do {
#define END } while (false)

// These preprocessor directives quietly disable GNU extensions
// if they are unsupported.
#if !defined __GNUC__ || __GNUC__ < 4
#ifndef __attribute__
#define __attribute__(_)
#endif
#endif

// The preprocessor directive `static_assert(p, s)`
// imitates the standard library function with the same name
// if it is not available.
// Due to technical limitations each `static_assert` must be on its own line
// to avoid naming conflicts.
#ifndef static_assert
#define _static_assert_line(p, n) __attribute__ ((__unused__)) \
  static int const _static_assert_##n[(p) ? 1 : -1]
#define _static_assert(p, n) _static_assert_line((p), n)
#define static_assert(p, _) _static_assert((p), __LINE__)
#endif

// The preprocessor directive `dynamic_assert(p, s)`
// is equivalent to `assert(p)` and
// exists just for the sake of consistency with `static_assert`.
#ifndef dynamic_assert
#define dynamic_assert(p, _) assert(p)
#endif

// These preprocessor directives ensure
// that exactly one of `NDEBUG` or `DEBUG` is defined.
#ifdef DEBUG
#ifdef NDEBUG
static_assert(false, "Contradictory debug directives");
#else
#endif
#else
#ifndef NDEBUG
#define NDEBUG
#endif
#endif

// The preprocessor directive `MIN(x, y)`
// expands to the lesser of `x` and `y`.
#define MIN(x, y) ((x) < (y) ? (x) : (y))

// The preprocessor directive `MAX(x, y)`
// expands to the greater of `x` and `y`.
#define MAX(x, y) ((x) > (y) ? (x) : (y))

#endif
