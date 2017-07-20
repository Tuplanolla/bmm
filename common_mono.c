#include <stddef.h>

extern inline void bmm_hsort(size_t,
    int (*)(size_t, size_t, void *),
    void (*)(size_t, size_t, void *), void *);
