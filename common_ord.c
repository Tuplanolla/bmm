#include <stdbool.h>
#include <stddef.h>

#include "ext.h"

extern inline int type(bmm_cmp, A)(A, A);

extern inline A type(bmm_clamp, A)(A, A, A);

extern inline A type(bmm_min, A)(A, A);

extern inline A type(bmm_max, A)(A, A);
