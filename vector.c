#include "ext.h"
#include "vector.h"

#define A double
#define D 2
#include "vector_fp.c"
#undef D
#undef A

#define A double
#define D 3
#include "vector_fp.c"
#undef D
#undef A
