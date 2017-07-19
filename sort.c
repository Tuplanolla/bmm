#include <stddef.h>

#include "sort.h"

extern inline void bmm_sort_heapsort(size_t,
    int (*)(size_t, size_t, void *),
    void (*)(size_t, size_t, void *), void *);
