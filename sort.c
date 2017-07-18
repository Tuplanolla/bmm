#include <stddef.h>

#include "sort.h"

extern inline void hsort(size_t,
    int (*)(size_t, size_t, void *),
    void (*)(size_t, size_t, void *), void *);
