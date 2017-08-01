#include "alias.h"
#include "ext.h"
#include "wrap.h"

extern inline float type(fabs, float)(float);

extern inline double type(fabs, double)(double);

extern inline long double type(fabs, long_double)(long double);

extern inline float type(trunc, float)(float);

extern inline double type(trunc, double)(double);

extern inline long double type(trunc, long_double)(long double);

extern inline float type(fmod, float)(float, float);

extern inline double type(fmod, double)(double, double);

extern inline long double type(fmod, long_double)(long double,
    long double);

extern inline float type(pow, float)(float, float);

extern inline double type(pow, double)(double, double);

extern inline long double type(pow, long_double)(long double,
    long double);
