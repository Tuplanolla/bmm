#include <stddef.h>

#include "ext.h"

extern inline $(bmm_quotrem_t, A) $(bmm_quotrem, A)(A, A);

extern inline A $(bmm_abs, A)(A);

extern inline A $(bmm_wrap, A)(A, A, A);

extern inline A $(bmm_uwrap, A)(A, A);

extern inline A $(bmm_uclamp, A)(A, A);

extern inline A $(bmm_uinc, A)(A, A);

extern inline A $(bmm_inc, A)(A, A, A);

extern inline A $(bmm_udec, A)(A, A);

extern inline A $(bmm_dec, A)(A, A, A);

extern inline A $(bmm_fact, A)(A);

extern inline A $(bmm_multfact, A)(A, A);

extern inline A $(bmm_tri, A)(A);

extern inline A $(bmm_tamean2, A)(A, A);

extern inline A $(bmm_famean2, A)(A, A);

extern inline A $(bmm_flog, A)(A, A);

extern inline A $(bmm_clog, A)(A, A);

extern inline A $(bmm_firt, A)(A, A);

extern inline A $(bmm_cirt, A)(A, A);

extern inline void $(bmm_hc, A)(A *, A,
    size_t, A);

extern inline A $(bmm_unhc, A)(A const *,
    size_t, A);

extern inline void $(bmm_hcd, A)(A *restrict, A,
    size_t, A const *restrict);

extern inline A $(bmm_unhcd, A)(A const *restrict,
    size_t, A const *restrict);
