#include <stdbool.h>
#include <stddef.h>

#include "neigh.h"

extern inline bool bmm_neigh_qijp(size_t*, size_t, size_t);

extern inline bool bmm_neigh_qij(size_t* restrict, size_t const* restrict,
    size_t, size_t, size_t const* restrict);

extern inline bool bmm_neigh_qijcp(size_t* restrict, size_t const* restrict,
    size_t, size_t, size_t const* restrict, bool const*);

extern inline size_t bmm_neigh_np(size_t);

extern inline size_t bmm_neigh_npr(size_t);

extern inline size_t bmm_neigh_nplh(size_t);

extern inline size_t bmm_neigh_nplhr(size_t);

extern inline size_t bmm_neigh_npuh(size_t);

extern inline size_t bmm_neigh_npuhr(size_t);

extern inline size_t bmm_neigh_n(size_t* restrict,
    size_t const* restrict, size_t, size_t const* restrict);

extern inline size_t bmm_neigh_nr(size_t* restrict,
    size_t const* restrict, size_t, size_t const* restrict);

extern inline size_t bmm_neigh_nlh(size_t* restrict,
    size_t const* restrict, size_t, size_t const* restrict);

extern inline size_t bmm_neigh_nlhr(size_t* restrict,
    size_t const* restrict, size_t, size_t const* restrict);

extern inline size_t bmm_neigh_nuh(size_t* restrict,
    size_t const* restrict, size_t, size_t const* restrict);

extern inline size_t bmm_neigh_nuhr(size_t* restrict,
    size_t const* restrict, size_t, size_t const* restrict);

extern inline size_t bmm_neigh_ncp(size_t* restrict,
    size_t const* restrict, size_t, size_t const* restrict, bool const*);

extern inline size_t bmm_neigh_ncpr(size_t* restrict,
    size_t const* restrict, size_t, size_t const* restrict, bool const*);

extern inline size_t bmm_neigh_ncplh(size_t* restrict,
    size_t const* restrict, size_t, size_t const* restrict, bool const*);

extern inline size_t bmm_neigh_ncplhr(size_t* restrict,
    size_t const* restrict, size_t, size_t const* restrict, bool const*);

extern inline size_t bmm_neigh_ncpuh(size_t* restrict,
    size_t const* restrict, size_t, size_t const* restrict, bool const*);

extern inline size_t bmm_neigh_ncpuhr(size_t* restrict,
    size_t const* restrict, size_t, size_t const* restrict, bool const*);

extern inline void bmm_neigh_ijp(size_t* restrict,
    size_t const* restrict, size_t, size_t, size_t const* restrict);

extern inline void bmm_neigh_ijpr(size_t* restrict,
    size_t const* restrict, size_t, size_t, size_t const* restrict);

extern inline void bmm_neigh_ijplh(size_t* restrict,
    size_t const* restrict, size_t, size_t, size_t const* restrict);

extern inline void bmm_neigh_ijplhr(size_t* restrict,
    size_t const* restrict, size_t, size_t, size_t const* restrict);

extern inline void bmm_neigh_ijpuh(size_t* restrict,
    size_t const* restrict, size_t, size_t, size_t const* restrict);

extern inline void bmm_neigh_ijpuhr(size_t* restrict,
    size_t const* restrict, size_t, size_t, size_t const* restrict);

extern inline void bmm_neigh_ij(size_t* restrict,
    size_t const* restrict, size_t, size_t, size_t const* restrict);

extern inline void bmm_neigh_ijr(size_t* restrict,
    size_t const* restrict, size_t, size_t, size_t const* restrict);

extern inline void bmm_neigh_ijlh(size_t* restrict,
    size_t const* restrict, size_t, size_t, size_t const* restrict);

extern inline void bmm_neigh_ijlhr(size_t* restrict,
    size_t const* restrict, size_t, size_t, size_t const* restrict);

extern inline void bmm_neigh_ijuh(size_t* restrict,
    size_t const* restrict, size_t, size_t, size_t const* restrict);

extern inline void bmm_neigh_ijuhr(size_t* restrict,
    size_t const* restrict, size_t, size_t, size_t const* restrict);

extern inline void bmm_neigh_ijcp(size_t* restrict, size_t const* restrict,
    size_t, size_t, size_t const* restrict, bool const*);

extern inline void bmm_neigh_ijcpr(size_t* restrict, size_t const* restrict,
    size_t, size_t, size_t const* restrict, bool const*);

extern inline void bmm_neigh_ijcplh(size_t* restrict, size_t const* restrict,
    size_t, size_t, size_t const* restrict, bool const*);

extern inline void bmm_neigh_ijcplhr(size_t* restrict, size_t const* restrict,
    size_t, size_t, size_t const* restrict, bool const*);

extern inline void bmm_neigh_ijcpuh(size_t* restrict, size_t const* restrict,
    size_t, size_t, size_t const* restrict, bool const*);

extern inline void bmm_neigh_ijcpuhr(size_t* restrict, size_t const* restrict,
    size_t, size_t, size_t const* restrict, bool const*);

extern inline size_t bmm_neigh_qicp(size_t const* restrict,
    size_t, size_t, size_t const* restrict, bool const*);

extern inline size_t bmm_neigh_icp(size_t const* restrict,
    size_t, size_t, size_t const* restrict, bool const*, int);

extern inline size_t bmm_neigh_icpuhr(size_t const* restrict,
    size_t, size_t, size_t const* restrict, bool const*);
