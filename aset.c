#include <stdbool.h>
#include <stddef.h>

#include "aset.h"

extern inline bool bmm_aset_canins(size_t*, size_t);

extern inline void bmm_aset_preins(void* restrict,
    size_t* restrict, size_t, void const* restrict);

extern inline void bmm_aset_ins(void* restrict,
    size_t* restrict, size_t, void const* restrict);

extern inline bool bmm_aset_candel(size_t*, size_t);

extern inline void bmm_aset_redel(void* restrict,
    size_t* restrict, size_t, size_t);

extern inline void bmm_aset_del(void* restrict,
    size_t* restrict, size_t, size_t);
