#ifndef BMM_SDL_H
/// Interactive viewer based on SDL and OpenGL.
#define BMM_SDL_H

#include <SDL2/SDL.h>
#include <limits.h>
#include <stdbool.h>
#include <stddef.h>
#include <sys/time.h>

#include "dem.h"
#include "ext.h"

/// The call `bmm_sdl_t_to_timeval(tp, t)`
/// sets the time structure `tp` to approximately `t` ticks.
__attribute__ ((__nonnull__))
inline void bmm_sdl_t_to_timeval(struct timeval* const tp, Uint32 const t) {
  tp->tv_sec = t / 1000;
  tp->tv_usec = t % 1000 * 1000;
}

/// The call `bmm_sdl_t_from_timeval(tp)`
/// returns the approximate number of ticks in the time structure `tp`.
__attribute__ ((__nonnull__, __pure__))
inline Uint32 bmm_sdl_t_from_timeval(struct timeval const* const tp) {
  return (Uint32) (tp->tv_sec * 1000 + tp->tv_usec / 1000);
}

/// The call `bmm_sdl_trem(tnow, tnext)`
/// returns the number of ticks from `tnow` to `tnext`.
/// The minimum image convention is applied to correct wrapping.
__attribute__ ((__const__, __pure__))
inline Uint32 bmm_sdl_trem(Uint32 const tnow, Uint32 const tnext) {
  if (tnow < tnext)
    return tnext - tnow;
  else {
    Uint32 const tdiff = tnow - tnext;
    return tdiff < UINT32_MAX / 2 ? 0 : UINT32_MAX - tdiff;
  }
}

struct bmm_sdl_opts {
  unsigned int width;
  unsigned int height;
  unsigned int fps;
  unsigned int ms;
  double zoomfac;
};

struct bmm_sdl {
  struct bmm_sdl_opts opts;
  int width;
  int height;
  double qaspect;
  double qzoom;
  double rorigin[2];
  unsigned int fps;
  bool stale;
  bool active;
  bool _pad0[2];
  size_t itarget;
  struct bmm_dem dem;
};

__attribute__ ((__nonnull__))
void bmm_sdl_opts_def(struct bmm_sdl_opts*);

__attribute__ ((__nonnull__))
void bmm_sdl_def(struct bmm_sdl*, struct bmm_sdl_opts const*);

__attribute__ ((__nonnull__))
bool bmm_sdl_run(struct bmm_sdl*);

__attribute__ ((__nonnull__))
bool bmm_sdl_run_with(struct bmm_sdl_opts const*);

#endif
