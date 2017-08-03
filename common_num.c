#include <stddef.h>

#include "ext.h"

extern inline int type(bmm_sgn, A)(A);

extern inline A type(bmm_abs, A)(A);

extern inline A type(bmm_power, A)(A, size_t);

extern inline A type(bmm_sum, A)(A const *, size_t);

extern inline A type(bmm_prod, A)(A const *, size_t);
