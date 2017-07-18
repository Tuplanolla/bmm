#include <stdbool.h>
#include <stddef.h>

#include "neigh.h"

extern inline bool bmm_neigh_qijcpij(size_t *restrict, size_t const *restrict,
    size_t, size_t, size_t const *restrict, bool const *);

extern inline bool bmm_neigh_qijcpi(size_t *restrict, size_t,
    size_t, size_t, size_t const *restrict, bool const *);

extern inline size_t bmm_neigh_qicpij(size_t const *restrict,
    size_t, size_t, size_t const *restrict, bool const *);

extern inline size_t bmm_neigh_qicpi(size_t,
    size_t, size_t, size_t const *, bool const *);

extern inline size_t bmm_neigh_ncpij(size_t const *restrict,
    size_t, size_t const *restrict, bool const *, int);

extern inline size_t bmm_neigh_ncpi(size_t,
    size_t, size_t const *, bool const *, int);

extern inline void bmm_neigh_ijcpij(size_t *restrict, size_t const *restrict,
    size_t, size_t, size_t const *restrict, bool const *, int);

extern inline void bmm_neigh_ijcpi(size_t *restrict, size_t,
    size_t, size_t, size_t const *restrict, bool const *, int);

extern inline size_t bmm_neigh_icpij(size_t const *restrict,
    size_t, size_t, size_t const *restrict, bool const *, int);

extern inline size_t bmm_neigh_icpi(size_t,
    size_t, size_t, size_t const *, bool const *, int);
