#include <stddef.h>

#include "ext.h"

extern inline int type(bmm_sgn, A)(A);

extern inline A type(bmm_min, A)(A, A);

extern inline A type(bmm_max, A)(A, A);

extern inline A type(bmm_abs, A)(A);

extern inline A type(bmm_pow, A)(A, size_t);
