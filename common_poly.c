#include <stddef.h>

#include "ext.h"

extern inline void $(bmm_swap, A)(A *restrict, A *restrict);

extern inline void $(bmm_map, A)(size_t, void (*)(size_t));

extern inline void $(bmm_map_cls, A)(size_t,
    void (*)(size_t, void *), void *);

extern inline A $(bmm_foldl, A)(size_t,
    A (*)(size_t, A), A);

extern inline A $(bmm_foldl_cls, A)(size_t,
    A (*)(size_t, A, void const *), A, void const *);

extern inline A $(bmm_foldr, A)(size_t,
    A (*)(size_t, A), A);

extern inline A $(bmm_foldr_cls, A)(size_t,
    A (*)(size_t, A, void const *), A, void const *);
