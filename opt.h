#ifndef BMM_OPT_H
/// Restricted command line option parsing.
#define BMM_OPT_H

#include <stdbool.h>
#include <stddef.h>

#include "ext.h"

/// The call `bmm_opt_parse(args, narg, f, ptr)`
/// parses the command line argument strings `args` of length `narg`.
/// When an argument pair in the form `--key value` is encountered,
/// the call `f("key", "value", ptr)` is made.
/// This allows parsing ordered dictionaries.
__attribute__ ((__nonnull__ (1, 3)))
bool bmm_opt_parse(char const* const*, size_t,
    bool (*)(char const*, char const*, void*), void*);

#endif
