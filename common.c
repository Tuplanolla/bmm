#include <stddef.h>

#include "alias.h"
#include "common.h"
#include "ext.h"

#include "common_mono.c"

#define A signed_char
#include "common_poly.c"
#undef A
#define A unsigned_char
#include "common_poly.c"
#undef A
#define A int
#include "common_poly.c"
#undef A
#define A double
#include "common_poly.c"
#undef A
#define A size_t
#include "common_poly.c"
#undef A

#define A signed_char
#include "common_ord.c"
#undef A
#define A unsigned_char
#include "common_ord.c"
#undef A
#define A int
#include "common_ord.c"
#undef A
#define A double
#include "common_ord.c"
#undef A
#define A size_t
#include "common_ord.c"
#undef A

#define A signed_char
#include "common_num.c"
#undef A
#define A unsigned_char
#include "common_num.c"
#undef A
#define A int
#include "common_num.c"
#undef A
#define A double
#include "common_num.c"
#undef A
#define A size_t
#include "common_num.c"
#undef A

#define A signed_char
#include "common_sint.c"
#undef A
#define A int
#include "common_sint.c"
#undef A

#define A unsigned_char
#include "common_uint.c"
#undef A
#define A size_t
#include "common_uint.c"
#undef A

#define A signed_char
#include "common_int.c"
#undef A
#define A unsigned_char
#include "common_int.c"
#undef A
#define A int
#include "common_int.c"
#undef A
#define A size_t
#include "common_int.c"
#undef A

#define A double
#include "common_fp.c"
#undef A
