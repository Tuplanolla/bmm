/// Synchronous signal handling.
#ifndef BMM_SIG_H
#define BMM_SIG_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#define BMM_SIG_MIN (SIG_ATOMIC_MIN + 2)
#define BMM_SIG_MAX (SIG_ATOMIC_MAX - 1)

/// The call `bmm_sig_handler(signum)`
/// handles the signal with the number `signum`.
/// This should be registered as a signal handler with `bmm_sig_register`.
void bmm_sig_handler(int);

/// The call `bmm_sig_register(sigs, n)` tries to register `bmm_sig_handler`
/// as the signal handler of each signal in the array `sigs` of length `n`.
/// If the operation is successful, `SIZE_MAX` is returned.
/// Otherwise the registering is stopped at the first failure and
/// the index of the first signal that failed to be registered is returned.
size_t bmm_sig_register(int const*, size_t);

/// The call `bmm_sig_unregister(sigs, n)` tries to restore the default
/// signal handler for each signal in the array `sigs` of length `n`.
/// If the operation is successful, `SIZE_MAX` is returned.
/// Otherwise the registering is stopped at the first failure and
/// the index of the first signal that failed to be unregistered is returned.
size_t bmm_sig_unregister(int const*, size_t);

/// The call `bmm_sig_unset()` returns `true` if a signal has not been caught.
/// Otherwise `false` is returned.
bool bmm_sig_unset(void);

/// The call `bmm_sig_set()` returns `true` if a signal has been caught.
/// Otherwise `false` is returned.
bool bmm_sig_set(void);

/// The call `bmm_sig_under()` returns `true`
/// if a signal below the representable minimum (negative) has been caught.
/// Otherwise `false` is returned.
bool bmm_sig_under(void);

/// The call `bmm_sig_normal(ptr)` returns `true`
/// if a normal signal (with a small enough magnitude) has been caught and
/// copies its value into `ptr` if `ptr` is not `NULL`.
/// Otherwise `false` is returned and no copying is done.
bool bmm_sig_normal(int*);

/// The call `bmm_sig_null()` returns `true`
/// if a null signal (zero) has been caught.
/// Otherwise `false` is returned.
bool bmm_sig_null(void);

/// The call `bmm_sig_over()` returns `true`
/// if a signal above the representable maximum (positive) has been caught.
/// Otherwise `false` is returned.
bool bmm_sig_over(void);

/// The call `bmm_sig_forget()` forgets about any previously caught signals.
/// Calling `bmm_sig_forget` before using any other signal handling procedures
/// from this compilation unit is not necessary.
void bmm_sig_forget(void);

/// The call `bmm_sig_use(ptr)` is equivalent
/// to performing `bmm_sig_normal(ptr)` and `bmm_sig_forget()` atomically.
bool bmm_sig_use(int*);

#endif
