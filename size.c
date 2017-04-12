#include <stdbool.h>
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

extern inline size_t bmm_size_firt(size_t, size_t);

extern inline size_t bmm_size_cirt(size_t, size_t);

extern inline size_t bmm_size_filog(size_t, size_t);

extern inline size_t bmm_size_cilog(size_t, size_t);

extern inline size_t bmm_size_uclamp(size_t, size_t);

extern inline size_t bmm_size_uwrap(size_t, size_t);

extern inline size_t bmm_size_inc(size_t, size_t);

extern inline size_t bmm_size_dec(size_t, size_t);

extern inline size_t bmm_size_sum(size_t const*, size_t);

extern inline size_t bmm_size_prod(size_t const*, size_t);

extern inline bool bmm_size_to_buffer(size_t*, unsigned char*, size_t,
    size_t, enum bmm_size_format);

extern inline bool bmm_size_from_buffer(size_t*, unsigned char const*, size_t,
    enum bmm_size_format);
