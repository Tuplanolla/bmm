#ifndef BMM_CLOPTS_H
#define BMM_CLOPTS_H

#include "exts.h"
#include <stdbool.h>
#include <stddef.h>

// The call `bmm_clopts(args, narg, f, p)`
// parses the command line argument strings `args` of length `narg`.
// When an argument pair in the form `--key value` is encountered,
// the call `f("key", "value", p)` is made.
__attribute__ ((__nonnull__ (1, 3)))
bool bmm_clopts(char const* const*, size_t,
    bool (*)(char const*, char const*, void*), void*);

#endif
