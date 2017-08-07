/// Vector arithemtic.

#ifndef BMM_VECTOR_H
#define BMM_VECTOR_H

#include <stddef.h>

#include "ext.h"
#include "scalar.h"

#define A double
#define D 2
#include "vector_fp.h"
#undef D
#undef A

#define A double
#define D 3
#include "vector_fp.h"
#undef D
#undef A

#endif
