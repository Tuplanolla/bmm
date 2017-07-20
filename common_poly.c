#include <stddef.h>

#include "ext.h"

extern inline void inst(bmm_swap, A)(A *restrict, A *restrict);

extern inline void inst(bmm_map, A)(size_t, void (*)(size_t));

extern inline void inst(bmm_map_cls, A, B)(B,
    void (*)(B, void *), void *);

extern inline A inst(bmm_foldl, A)(size_t,
    A (*)(size_t, A), A);

extern inline A inst(bmm_foldl_cls, A)(size_t,
    A (*)(size_t, A, void *), A, void *);

extern inline A inst(bmm_foldr, A)(size_t,
    A (*)(size_t, A), A);

extern inline A inst(bmm_foldr_cls, A)(size_t,
    A (*)(size_t, A, void *), A, void *);
