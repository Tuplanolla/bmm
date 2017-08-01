#include "alias.h"
#include "ext.h"
#include "wrap.h"

extern inline float type(bmm_trunc, float)(float);

extern inline double type(bmm_trunc, double)(double);

extern inline long double type(bmm_trunc, long_double)(long double);

extern inline float type(bmm_fmod, float)(float, float);

extern inline double type(bmm_fmod, double)(double, double);

extern inline long double type(bmm_fmod, long_double)(long double,
    long double);

extern inline float type(bmm_pow, float)(float, float);

extern inline double type(bmm_pow, double)(double, double);

extern inline long double type(bmm_pow, long_double)(long double,
    long double);
