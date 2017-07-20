#include <stddef.h>

#include "common.h"
#include "ext.h"

#include "common_mono.c"

#define B size_t
#define A int
#include "common_poly.c"
#undef A

#define A double
#include "common_poly.c"
#undef A

#define A size_t
#include "common_poly.c"
#undef A

#define A int
#include "common_ord.c"
#undef A
