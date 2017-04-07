/// Approximate time conversions using seconds.
#ifndef BMM_SEC_H
#define BMM_SEC_H

#include <math.h>
#include <sys/time.h>
#include <time.h>

#include "ext.h"

/// The call `bmm_sec_to_timeval(tp, s)`
/// sets the time structure `tp` to approximately `s` seconds.
__attribute__ ((__nonnull__))
inline void bmm_sec_to_timeval(struct timeval* const tp, double const s) {
  long int const sec = (long int) s;

  tp->tv_sec = sec;
  tp->tv_usec = (long int) ((s - (double) sec) * 1.0e+6);
}

/// The call `bmm_sec_from_timeval(tp)`
/// returns the approximate number of seconds in the time structure `tp`.
__attribute__ ((__nonnull__, __pure__))
inline double bmm_sec_from_timeval(struct timeval const* const tp) {
  return (double) tp->tv_sec + (double) tp->tv_usec / 1.0e+6;
}

/// The call `bmm_sec_to_timespec(tp, s)`
/// sets the time structure `tp` to approximately `s` seconds.
__attribute__ ((__nonnull__))
inline void bmm_sec_to_timespec(struct timespec* const tp, double const s) {
  long int const sec = (long int) s;

  tp->tv_sec = (time_t) sec;
  tp->tv_nsec = (long int) ((s - (double) sec) * 1.0e+9);
}

/// The call `bmm_sec_from_timespec(tp)`
/// returns the approximate number of seconds in the time structure `tp`.
__attribute__ ((__nonnull__, __pure__))
inline double bmm_sec_from_timespec(struct timespec const* const tp) {
  return (double) tp->tv_sec + (double) tp->tv_nsec / 1.0e+9;
}

/// The call `bmm_sec_now()`
/// returns the approximate monotonic time in seconds right now if possible.
/// Otherwise `NAN` is returned.
inline double bmm_sec_now(void) {
  struct timespec tp;
  return clock_gettime(CLOCK_MONOTONIC, &tp) == -1 ?
    (double) NAN : bmm_sec_from_timespec(&tp);
}

#endif
