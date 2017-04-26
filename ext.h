#ifndef BMM_EXT_H
/// Language extensions.
///
/// Note that this header defines symbols beyond its own namespace,
/// so it header should not be included before standard library headers.
#define BMM_EXT_H

#include <assert.h>
#include <stdbool.h>

/// These preprocessor directives quietly disable GNU extensions
/// if they are unsupported.
#if !defined __GNUC__ || __GNUC__ < 4

#ifndef __attribute__
#define __attribute__(_)
#endif

#endif

/// The preprocessor directive `thread_local`
/// expands to the optimal thread-local storage specifier.
#ifdef _GNU_SOURCE
#define thread_local __thread
#else
#define thread_local _Thread_local
#endif

/// The preprocessor directive `ever` makes it possible to write `for ever`.
#ifndef ever
#define ever (;;)
#endif

/// The preprocessor directives `begin` and `end` allow converting
/// multiple statements into a single statement
/// with local scope and no redundant semicolons.
/// They must always appear together.
#define begin do {
#define end } while (false)

/// The preprocessor directive `nmembof(x)`
/// expands to the number of elements in the array `x`.
#define nmembof(x) (sizeof (x) / sizeof *(x))

/// The preprocessor directive `static_assert(p, s)`
/// imitates the standard library function with the same name
/// if it is not available.
/// Due to technical limitations each `static_assert` must be on its own line
/// to avoid naming conflicts.
#ifndef static_assert
#define static_assert_line_(p, n) __attribute__ ((__unused__)) \
  static int const _static_assert_##n[(p) ? 1 : -1]
#define static_assert_(p, n) static_assert_line_((p), n)
#define static_assert(p, _) static_assert_((p), __LINE__)
#endif

/// The preprocessor directive `dynamic_assert(p, s)`
/// is equivalent to `assert(p)` with the comment `s`.
/// It exists for the sake of consistency with `static_assert`.
#ifndef dynamic_assert
#define dynamic_assert(p, _) assert(p)
#endif

/// These preprocessor directives ensure
/// that exactly one of `NDEBUG` or `DEBUG` is defined.
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

#endif
