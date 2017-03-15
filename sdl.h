// Interactive viewer based on SDL and OpenGL.
#ifndef BMM_SDL_H
#define BMM_SDL_H

#include "ext.h"
#include <SDL/SDL.h>
#include <stdbool.h>
#include <stddef.h>
#include <sys/time.h>

struct bmm_sdl_opts {
  size_t width;
  size_t height;
  unsigned int fps;
  unsigned int msaa;
};

// The call `bmm_sdl_ticks_to_timeval(tp, t)` sets
// the time structure `tp` to approximately `t` ticks.
__attribute__ ((__nonnull__))
inline void bmm_sdl_ticks_to_timeval(struct timeval* const tp,
    Uint32 const t) {
  tp->tv_sec = t / 1000;
  tp->tv_usec = t % 1000 * 1000;
}

// The call `bmm_sdl_ticks_from_timeval(tp)` returns
// the approximate number of ticks in the time structure `tp`.
__attribute__ ((__nonnull__, __pure__))
inline Uint32 bmm_sdl_ticks_from_timeval(struct timeval const* const tp) {
  return (Uint32) (tp->tv_sec * 1000 + tp->tv_usec / 1000);
}

__attribute__ ((__nonnull__))
void bmm_sdl_defopts(struct bmm_sdl_opts*);

__attribute__ ((__nonnull__))
bool bmm_sdl_run(struct bmm_sdl_opts const*);

#endif
