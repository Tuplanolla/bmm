#include <stddef.h>

#include "ext.h"

extern inline type(bmm_quot_t, A) type(bmm_quot, A)(A, A);

extern inline A type(bmm_abs, A)(A);

extern inline A type(bmm_wrap, A)(A, A, A);

extern inline A type(bmm_uwrap, A)(A, A);

extern inline A type(bmm_fact, A)(A);

extern inline A type(bmm_multfact, A)(A, A);

extern inline A type(bmm_tamean2, A)(A, A);

extern inline A type(bmm_famean2, A)(A, A);

extern inline A type(bmm_flog, A)(A, A);

extern inline A type(bmm_clog, A)(A, A);

extern inline void type(bmm_hc, A)(A *, A,
    size_t, A);

extern inline A type(bmm_unhc, A)(A const *,
    size_t, A);

extern inline void type(bmm_hcd, A)(A *restrict, A,
    size_t, A const *restrict);

extern inline A type(bmm_unhcd, A)(A const *restrict,
    size_t, A const *restrict);
