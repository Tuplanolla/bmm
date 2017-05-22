#include <stdbool.h>
#include <stddef.h>

#include "moore.h"

extern inline bool bmm_moore_qp(size_t*, size_t, size_t);

extern inline bool bmm_moore_q(size_t* restrict, size_t const* restrict,
    size_t, size_t, size_t const* restrict);

extern inline bool bmm_moore_qcp(size_t* restrict, size_t const* restrict,
    size_t, size_t, size_t const* restrict, bool const*);

extern inline size_t bmm_moore_np(size_t);

extern inline size_t bmm_moore_npr(size_t);

extern inline size_t bmm_moore_nplh(size_t);

extern inline size_t bmm_moore_nplhr(size_t);

extern inline size_t bmm_moore_npuh(size_t);

extern inline size_t bmm_moore_npuhr(size_t);

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

extern inline size_t bmm_moore_ncp(size_t* restrict,
    size_t const* restrict, size_t, size_t const* restrict, bool const*);

extern inline size_t bmm_moore_ncpr(size_t* restrict,
    size_t const* restrict, size_t, size_t const* restrict, bool const*);

extern inline size_t bmm_moore_ncplh(size_t* restrict,
    size_t const* restrict, size_t, size_t const* restrict, bool const*);

extern inline size_t bmm_moore_ncplhr(size_t* restrict,
    size_t const* restrict, size_t, size_t const* restrict, bool const*);

extern inline size_t bmm_moore_ncpuh(size_t* restrict,
    size_t const* restrict, size_t, size_t const* restrict, bool const*);

extern inline size_t bmm_moore_ncpuhr(size_t* restrict,
    size_t const* restrict, size_t, size_t const* restrict, bool const*);

extern inline void bmm_moore_ijp(size_t* restrict,
    size_t const* restrict, size_t, size_t, size_t const* restrict);

extern inline void bmm_moore_ijpr(size_t* restrict,
    size_t const* restrict, size_t, size_t, size_t const* restrict);

extern inline void bmm_moore_ijplh(size_t* restrict,
    size_t const* restrict, size_t, size_t, size_t const* restrict);

extern inline void bmm_moore_ijplhr(size_t* restrict,
    size_t const* restrict, size_t, size_t, size_t const* restrict);

extern inline void bmm_moore_ijpuh(size_t* restrict,
    size_t const* restrict, size_t, size_t, size_t const* restrict);

extern inline void bmm_moore_ijpuhr(size_t* restrict,
    size_t const* restrict, size_t, size_t, size_t const* restrict);

extern inline void bmm_moore_ij(size_t* restrict,
    size_t const* restrict, size_t, size_t, size_t const* restrict);

extern inline void bmm_moore_ijr(size_t* restrict,
    size_t const* restrict, size_t, size_t, size_t const* restrict);

extern inline void bmm_moore_ijlh(size_t* restrict,
    size_t const* restrict, size_t, size_t, size_t const* restrict);

extern inline void bmm_moore_ijlhr(size_t* restrict,
    size_t const* restrict, size_t, size_t, size_t const* restrict);

extern inline void bmm_moore_ijuh(size_t* restrict,
    size_t const* restrict, size_t, size_t, size_t const* restrict);

extern inline void bmm_moore_ijuhr(size_t* restrict,
    size_t const* restrict, size_t, size_t, size_t const* restrict);

extern inline void bmm_moore_ijcp(size_t* restrict, size_t const* restrict,
    size_t, size_t, size_t const* restrict, bool const*);

extern inline void bmm_moore_ijcpr(size_t* restrict, size_t const* restrict,
    size_t, size_t, size_t const* restrict, bool const*);

extern inline void bmm_moore_ijcplh(size_t* restrict, size_t const* restrict,
    size_t, size_t, size_t const* restrict, bool const*);

extern inline void bmm_moore_ijcplhr(size_t* restrict, size_t const* restrict,
    size_t, size_t, size_t const* restrict, bool const*);

extern inline void bmm_moore_ijcpuh(size_t* restrict, size_t const* restrict,
    size_t, size_t, size_t const* restrict, bool const*);

extern inline void bmm_moore_ijcpuhr(size_t* restrict, size_t const* restrict,
    size_t, size_t, size_t const* restrict, bool const*);
