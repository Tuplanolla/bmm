#include <stddef.h>

#include "ext.h"

extern inline inst(bmm_div_t, A) inst(bmm_div, A)(A, A);

extern inline void inst(bmm_hc, A)(A *, A,
    size_t, A);

extern inline A inst(bmm_unhc, A)(A const *,
    size_t, A);

extern inline void inst(bmm_hcd, A)(A *restrict, A,
    size_t, A const *restrict);

extern inline A inst(bmm_unhcd, A)(A const *restrict,
    size_t, A const *restrict);
