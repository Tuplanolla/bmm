#include <stdbool.h>
#include <stddef.h>

#include "alias.h"
#include "ext.h"
#include "scalar.h"

#include "scalar_lift.c"

#define A signed_char
#include "scalar_int.c"
#undef A
#define A int
#include "scalar_int.c"
#undef A

#define A signed_char
#include "scalar_sint.c"
#undef A
#define A int
#include "scalar_sint.c"
#undef A

#define A unsigned_char
#include "scalar_uint.c"
#undef A
#define A size_t
#include "scalar_uint.c"
#undef A

#define A double
#include "scalar_fp.c"
#undef A
