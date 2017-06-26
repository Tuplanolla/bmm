#include <stdbool.h>
#include <stddef.h>

#include "neigh.h"

extern inline bool bmm_neigh_qijcp(size_t* restrict, size_t const* restrict,
    size_t, size_t, size_t const* restrict, bool const*);

extern inline size_t bmm_neigh_qicp(size_t const* restrict,
    size_t, size_t, size_t const* restrict, bool const*);

extern inline size_t bmm_neigh_qicplin(size_t,
    size_t, size_t, size_t const*, bool const*);

extern inline size_t bmm_neigh_ncp(size_t const* restrict,
    size_t, size_t const* restrict, bool const*, int);

extern inline void bmm_neigh_ijcp(size_t* restrict, size_t const* restrict,
    size_t, size_t, size_t const* restrict, bool const*, int);

extern inline size_t bmm_neigh_icp(size_t const* restrict,
    size_t, size_t, size_t const* restrict, bool const*, int);

extern inline size_t bmm_neigh_icplin(size_t,
    size_t, size_t, size_t const*, bool const*, int);
