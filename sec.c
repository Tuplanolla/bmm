#include <sys/time.h>
#include <time.h>

#include "sec.h"

extern inline void bmm_sec_to_timeval(struct timeval*, double);

extern inline double bmm_sec_from_timeval(struct timeval const*);

extern inline void bmm_sec_to_timespec(struct timespec*, double);

extern inline double bmm_sec_from_timespec(struct timespec const*);

extern inline double bmm_sec_now(void);
