#ifndef BMM_HACK_H
/// Ugly hacks.
#define BMM_HACK_H

#include <stdbool.h>
#include <stddef.h>

/// The call `bmm_hack_strerror_r(num, buf, n)`
/// is equivalent to `strerror_r(num, buf, n) == 0`
/// with traditional `errno` handling and
/// XSI compliance instead of GNU compliance.
bool bmm_hack_strerror_r(int, char *, size_t);

#endif
