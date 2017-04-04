// Messaging protocol.
#ifndef BMM_MSG_H
#define BMM_MSG_H

#include "ext.h"
#include "dem.h"
#include <stdbool.h>
#include <stdio.h>

struct bmm_msg_head {
  unsigned char flags;
  unsigned char type;
};

#define BMM_FBIT_INTLE 7
#define BMM_FBIT_FPLE 5
#define BMM_FBIT_FLUSH 4
#define BMM_FBIT_BODY 3
#define BMM_FBIT_PREFIX 2

enum bmm_msg {
  BMM_MSG_NOP = 0,
  BMM_MSG_NPART = 42,
  BMM_MSG_NSTEP = 43,
  BMM_MSG_PARTS = 44,
  BMM_MSG_EKINE = 45,
  BMM_MSG_ESMOM = 46,
  BMM_MSG_EVMOM = 47,
  BMM_MSG_NEIGH = 48,
  BMM_MSG_URGH = 80
};

static_assert(BMM_MSG_URGH < BMM_MSG_MAX, "Too many messages");

// TODO Wrangle prototypes.

__attribute__ ((__nonnull__ (2)))
bool bmm_msg_preread(size_t* const ptr, struct bmm_msg_head const* const head,
    size_t const size);

__attribute__ ((__nonnull__))
bool bmm_msg_prewrite(struct bmm_msg_head const* const bad_head,
    size_t const size);

__attribute__ ((__nonnull__))
void bmm_defhead(struct bmm_msg_head*);

__attribute__ ((__deprecated__, __nonnull__))
bool bmm_msg_get(struct bmm_msg_head*, struct bmm_dem*);

__attribute__ ((__deprecated__, __nonnull__))
void bmm_msg_put(struct bmm_msg_head const*, struct bmm_dem const*);

#endif
