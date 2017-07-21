#include <stddef.h>

#include "ext.h"

extern inline A type(bmm_wrap, A)(A, A, A);

extern inline void type(bmm_hc, A)(A *, A,
    size_t, A);

extern inline A type(bmm_unhc, A)(A const *,
    size_t, A);

extern inline void type(bmm_hcd, A)(A *restrict, A,
    size_t, A const *restrict);

extern inline A type(bmm_unhcd, A)(A const *restrict,
    size_t, A const *restrict);
