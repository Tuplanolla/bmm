#include <stdbool.h>
#include <stddef.h>

#include "astack.h"

extern inline bool bmm_astack_canins(size_t*, size_t);

extern inline void bmm_astack_preins(void* restrict,
    size_t* restrict, size_t, void const* restrict);

extern inline void bmm_astack_ins(void* restrict,
    size_t* restrict, size_t, void const* restrict);

extern inline bool bmm_astack_candel(size_t*);

extern inline void bmm_astack_del(void* restrict,
    size_t* restrict, size_t);
