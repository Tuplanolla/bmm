// Messaging protocol.
#ifndef BMM_MSG_H
#define BMM_MSG_H

#include "dem.h"
#include <stdbool.h>
#include <stdio.h>

struct bmm_msg_head {
  unsigned char flags;
  unsigned char type;
};

#define BMM_FBIT_INTLE 7
#define BMM_FBIT_FPLE 6
#define BMM_FBIT_FPFMT 5
#define BMM_FBIT_FLUSH 4
#define BMM_FBIT_VARZ 3
#define BMM_FBIT_ZPRE 2

enum bmm_msg {
  BMM_MSG_NOP = 0x0,
  BMM_MSG_NPART = 0x42,
  BMM_MSG_NSTEP = 0x43,
  BMM_MSG_PARTS = 0x44,
  BMM_MSG_EKINE = 0x45,
  BMM_MSG_ESMOM = 0x46,
  BMM_MSG_EVMOM = 0x47,
  BMM_MSG_NEIGH = 0x48,
  BMM_MSG_URGH = 0x80
};

__attribute__ ((__nonnull__))
void bmm_defhead(struct bmm_msg_head*);

// Read a message and call `f` on the header.
// Return `false` in case of an error.
// Return `true` otherwise.
// With `EOF` return `true`, but do not call `f`.
// TODO Maybe change this?
__attribute__ ((__nonnull__ (1, 2, 3)))
bool bmm_msg_get(struct bmm_msg_head*, struct bmm_dem*,
    bool (*)(struct bmm_msg_head const*, void*), void*);

__attribute__ ((__nonnull__))
void bmm_msg_put(struct bmm_msg_head const*, struct bmm_dem const*);

#endif
