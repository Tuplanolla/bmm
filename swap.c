#include <stdbool.h>
#include <stddef.h>

#include "ext.h"
#include "swap.h"

#define T int
#include "swap_meta.c"
#undef T

#define T size_t
#include "swap_meta.c"
#undef T

#define T double
#include "swap_meta.c"
#undef T

#define T unsigned_int
#include "swap_meta.c"
#undef T
