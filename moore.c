#include <stdbool.h>
#include <stddef.h>

#include "moore.h"

extern inline bool bmm_moore_qp(void);

extern inline size_t bmm_moore_np(size_t);

extern inline size_t bmm_moore_npr(size_t);

extern inline size_t bmm_moore_nph(size_t);

extern inline size_t bmm_moore_nphr(size_t);

extern inline bool bmm_moore_q(size_t* restrict,
    size_t const* restrict, size_t, size_t, size_t const* restrict);

extern inline size_t bmm_moore_n(size_t* restrict,
    size_t const* restrict, size_t, size_t const* restrict);

extern inline size_t bmm_moore_nr(size_t* restrict,
    size_t const* restrict, size_t, size_t const* restrict);

extern inline size_t bmm_moore_nlh(size_t* restrict,
    size_t const* restrict, size_t, size_t const* restrict);

extern inline size_t bmm_moore_nlhr(size_t* restrict,
    size_t const* restrict, size_t, size_t const* restrict);

extern inline size_t bmm_moore_nuh(size_t* restrict,
    size_t const* restrict, size_t, size_t const* restrict);

extern inline size_t bmm_moore_nuhr(size_t* restrict,
    size_t const* restrict, size_t, size_t const* restrict);

extern inline bool bmm_moore_qcp(size_t* restrict,
    size_t const* restrict, size_t, size_t, size_t const* restrict);

extern inline size_t bmm_moore_ncp(size_t* restrict,
    size_t const* restrict, size_t, size_t const* restrict);

extern inline size_t bmm_moore_ncpr(size_t* restrict,
    size_t const* restrict, size_t, size_t const* restrict);

extern inline size_t bmm_moore_ncplh(size_t* restrict,
    size_t const* restrict, size_t, size_t const* restrict);

extern inline size_t bmm_moore_ncplhr(size_t* restrict,
    size_t const* restrict, size_t, size_t const* restrict);

extern inline size_t bmm_moore_ncpuh(size_t* restrict,
    size_t const* restrict, size_t, size_t const* restrict);

extern inline size_t bmm_moore_ncpuhr(size_t* restrict,
    size_t const* restrict, size_t, size_t const* restrict);

extern inline void bmm_moore_ij(size_t* restrict,
    size_t const* restrict, size_t, size_t);

extern inline void bmm_moore_ijr(size_t* restrict,
    size_t const* restrict, size_t, size_t);

extern inline void bmm_moore_ijh(size_t* restrict,
    size_t const* restrict, size_t, size_t);

extern inline void bmm_moore_ijhr(size_t* restrict,
    size_t const* restrict, size_t, size_t);

extern inline size_t bmm_moore_i(size_t const*,
    size_t, size_t);
