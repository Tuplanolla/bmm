#ifndef BMM_COMMON_H
/// Common operations.
#define BMM_COMMON_H

#include <stddef.h>

#include "ext.h"

#include "common_mono.h"

#define A int
#include "common_poly.h"
#undef A

#define A double
#include "common_poly.h"
#undef A

#define A size_t
#include "common_poly.h"
#undef A

#define A int
#include "common_ord.h"
#undef A

#endif
