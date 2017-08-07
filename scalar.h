/// Scalar arithemtic.

#ifndef BMM_SCALAR_H
#define BMM_SCALAR_H

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "alias.h"
#include "ext.h"

#include "scalar_lift.h"

#define A signed_char
#include "scalar_int.h"
#undef A
#define A int
#include "scalar_int.h"
#undef A

#define A signed_char
#include "scalar_sint.h"
#undef A
#define A int
#include "scalar_sint.h"
#undef A

#define A unsigned_char
#include "scalar_uint.h"
#undef A
#define A size_t
#include "scalar_uint.h"
#undef A

#define A double
#include "scalar_fp.h"
#undef A

#endif
