#include "sig.h"
#include <limits.h>
#include <signal.h>
#include <stdbool.h>
#include <stddef.h>

#define BMM_SIG_UNSET (BMM_SIG_MIN - 2)
#define BMM_SIG_UNDER (BMM_SIG_MIN - 1)
#define BMM_SIG_NULL 0
#define BMM_SIG_OVER (BMM_SIG_MAX + 1)

static volatile sig_atomic_t sig = BMM_SIG_UNSET;

static bool signum_under(int const signum) {
  return signum < BMM_SIG_MIN;
}

static bool signum_normal(int const signum) {
  return signum >= BMM_SIG_MIN && signum <= BMM_SIG_MAX;
}

static bool signum_null(int const signum) {
  return signum == 0;
}

static bool signum_over(int const signum) {
  return signum > BMM_SIG_MAX;
}

void bmm_sig_handler(int const signum) {
  sig = signum_normal(signum) ? (sig_atomic_t) signum :
    signum_null(signum) ? BMM_SIG_NULL :
    signum_over(signum) ? BMM_SIG_OVER :
    signum_under(signum) ? BMM_SIG_UNDER :
    BMM_SIG_UNSET;
}

// This is also possible to accomplish with `signal`
// if it provides BSD semantics instead of System V semantics,
// but why bother when `sigaction` exists?
size_t bmm_sig_register(int const* const sigs, size_t const n) {
  struct sigaction act;
  act.sa_handler = bmm_sig_handler;
  sigfillset(&act.sa_mask);
  act.sa_flags = SA_RESTART;

  for (size_t i = 0; i < n; ++i)
    if (sigaction(sigs[i], &act, NULL) == -1)
      return i;

  return SIZE_MAX;
}

size_t bmm_sig_unregister(int const* const sigs, size_t const n) {
  struct sigaction act;
  act.sa_handler = SIG_DFL;
  sigfillset(&act.sa_mask);
  act.sa_flags = 0;

  for (size_t i = 0; i < n; ++i)
    if (sigaction(sigs[i], &act, NULL) == -1)
      return i;

  return SIZE_MAX;
}

bool bmm_sig_unset(void) {
  return sig == BMM_SIG_UNSET;
}

bool bmm_sig_set(void) {
  return sig != BMM_SIG_UNSET;
}

bool bmm_sig_under(void) {
  return sig == BMM_SIG_MIN;
}

bool bmm_sig_normal(int* const ptr) {
  sig_atomic_t const n = sig;

  bool const p = n >= BMM_SIG_MIN && n <= BMM_SIG_MAX;

  if (p && ptr != NULL)
    *ptr = (int) n;

  return p;
}

bool bmm_sig_null(void) {
  return sig == BMM_SIG_NULL;
}

bool bmm_sig_over(void) {
  return sig == BMM_SIG_MAX;
}

void bmm_sig_forget(void) {
  sig = BMM_SIG_UNSET;
}

bool bmm_sig_use(int* const ptr) {
  sigset_t set;
  sigfillset(&set);

  sigset_t oldset;
  bool const p = sigprocmask(SIG_SETMASK, &set, &oldset) == -1;

  bool const q = bmm_sig_normal(ptr);
  if (q)
    bmm_sig_forget();

  if (!p)
    (void) sigprocmask(SIG_SETMASK, &oldset, NULL);

  return q;
}
