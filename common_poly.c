#include <stddef.h>

#include "ext.h"

extern inline void type(bmm_swap, A)(A *restrict, A *restrict);

extern inline void type(bmm_map, A)(size_t, void (*)(size_t));

extern inline void type(bmm_map_cls, A)(size_t,
    void (*)(size_t, void *), void *);

extern inline A type(bmm_foldl, A)(size_t,
    A (*)(size_t, A), A);

extern inline A type(bmm_foldl_cls, A)(size_t,
    A (*)(size_t, A, void *), A, void *);

extern inline A type(bmm_foldr, A)(size_t,
    A (*)(size_t, A), A);

extern inline A type(bmm_foldr_cls, A)(size_t,
    A (*)(size_t, A, void *), A, void *);
