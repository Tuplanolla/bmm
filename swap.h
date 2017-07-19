#ifndef BMM_SWAP_H
/// Swapping (test).
#define BMM_SWAP_H

#include <stdbool.h>
#include <stddef.h>

#include "ext.h"

#define T int
#include "swap_meta.h"
#undef T

#define T size_t
#include "swap_meta.h"
#undef T

#define T double
#include "swap_meta.h"
#undef T

// TODO Just to have these on record...
typedef signed char signed_char;
typedef unsigned char unsigned_char;
typedef short int short_int;
typedef unsigned short int unsigned_short_int;
typedef unsigned int unsigned_int;
typedef long int long_int;
typedef unsigned long int unsigned_long_int;
typedef long long int long_long_int;
typedef unsigned long long int unsigned_long_long_int;

#define T unsigned_int
#include "swap_meta.h"
#undef T

#endif
