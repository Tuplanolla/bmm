/// Common operations.

#ifndef BMM_COMMON_H
#define BMM_COMMON_H

#include <math.h>
#include <stddef.h>

#include "alias.h"
#include "ext.h"

#include "common_mono.h"

#define A signed_char
#include "common_poly.h"
#undef A
#define A unsigned_char
#include "common_poly.h"
#undef A
#define A int
#include "common_poly.h"
#undef A
#define A double
#include "common_poly.h"
#undef A
#define A size_t
#include "common_poly.h"
#undef A

#define A signed_char
#include "common_ord.h"
#undef A
#define A unsigned_char
#include "common_ord.h"
#undef A
#define A int
#include "common_ord.h"
#undef A
#define A double
#include "common_ord.h"
#undef A
#define A size_t
#include "common_ord.h"
#undef A

#define A signed_char
#include "common_num.h"
#undef A
#define A unsigned_char
#include "common_num.h"
#undef A
#define A int
#include "common_num.h"
#undef A
#define A double
#include "common_num.h"
#undef A
#define A size_t
#include "common_num.h"
#undef A

#define A signed_char
#include "common_int.h"
#undef A
#define A unsigned_char
#include "common_int.h"
#undef A
#define A int
#include "common_int.h"
#undef A
#define A size_t
#include "common_int.h"
#undef A

#define A signed_char
#include "common_sint.h"
#undef A
#define A int
#include "common_sint.h"
#undef A

#define A unsigned_char
#include "common_uint.h"
#undef A
#define A size_t
#include "common_uint.h"
#undef A

#define A double
#define TRUNCA trunc
#define FMODA fmod
#include "common_fp.h"
#undef FMODA
#undef TRUNCA
#undef A

#endif
