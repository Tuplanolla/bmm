#include <stddef.h>

#include "size.h"

extern inline bmm_size_div_t bmm_size_div(size_t, size_t);

extern inline int bmm_size_cmp(size_t, size_t);

extern inline size_t bmm_size_min(size_t, size_t);

extern inline size_t bmm_size_max(size_t, size_t);

extern inline size_t bmm_size_pow(size_t, size_t);

extern inline size_t bmm_size_identity(size_t);

extern inline size_t bmm_size_constant(size_t, size_t);

extern inline size_t bmm_size_zero(size_t);

extern inline size_t bmm_size_one(size_t);

extern inline size_t bmm_size_midpoint(size_t, size_t);

extern inline size_t bmm_size_sq(size_t);

extern inline size_t bmm_size_cb(size_t);

extern inline size_t bmm_size_firt(size_t, size_t);

extern inline size_t bmm_size_cirt(size_t, size_t);

extern inline size_t bmm_size_flog(size_t, size_t);

extern inline size_t bmm_size_clog(size_t, size_t);

extern inline size_t bmm_size_uclamp(size_t, size_t);

extern inline size_t bmm_size_wrap(size_t, size_t, size_t);

extern inline size_t bmm_size_uwrap(size_t, size_t);

extern inline size_t bmm_size_uinc(size_t, size_t);

extern inline size_t bmm_size_inc(size_t, size_t, size_t);

extern inline size_t bmm_size_udec(size_t, size_t);

extern inline size_t bmm_size_dec(size_t, size_t, size_t);

extern inline size_t bmm_size_fact(size_t);

extern inline size_t bmm_size_sum(size_t const*, size_t);

extern inline size_t bmm_size_prod(size_t const*, size_t);

extern inline size_t bmm_size_lfold(size_t (*)(size_t, size_t, void*),
    size_t const* restrict, size_t, size_t, void* restrict);

extern inline size_t bmm_size_rfold(size_t (*)(size_t, size_t, void*),
    size_t const* restrict, size_t, size_t, void* restrict);

extern inline void bmm_size_hc(size_t*, size_t, size_t, size_t);

extern inline size_t bmm_size_unhc(size_t const*, size_t, size_t);

extern inline void bmm_size_hcd(size_t* restrict,
    size_t, size_t, size_t const* restrict);

extern inline size_t bmm_size_unhcd(size_t const* restrict,
    size_t, size_t const* restrict);
