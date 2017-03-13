// Approximate time conversions from and to seconds.
#ifndef BMM_SEC_H
#define BMM_SEC_H

#include "ext.h"
#include <sys/time.h>
#include <time.h>

// The call `bmm_sec_to_timeval(tp, s)` sets
// the time structure `tp` to approximately `s` seconds.
__attribute__ ((__nonnull__, __pure__))
void bmm_sec_to_timeval(struct timeval*, double);

// The call `bmm_sec_from_timeval(tp)` returns
// the approximate number of seconds in the time structure `tp`.
__attribute__ ((__nonnull__, __pure__))
double bmm_sec_from_timeval(struct timeval const*);

// The call `bmm_sec_to_timespec(tp, s)` sets
// the time structure `tp` to approximately `s` seconds.
__attribute__ ((__nonnull__, __pure__))
void bmm_sec_to_timespec(struct timespec*, double);

// The call `bmm_sec_from_timespec(tp)` returns
// the approximate number of seconds in the time structure `tp`.
__attribute__ ((__nonnull__, __pure__))
double bmm_sec_from_timespec(struct timespec const*);

// The call `bmm_sec_now()` returns
// the approximate monotonic time in seconds right now if possible.
// Otherwise `NAN` is returned.
double bmm_sec_now(void);

#endif
