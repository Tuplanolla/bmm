// Messaging protocol.
#ifndef BMM_MSG_H
#define BMM_MSG_H

#include "dem.h"
#include <stdbool.h>
#include <stdio.h>

struct bmm_head {
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
  BMM_MSG_NDIM = 0x40,
  BMM_MSG_NBIN = 0x41,
  BMM_MSG_NPART = 0x42,
  BMM_MSG_NSTEP = 0x43,
  BMM_MSG_PARTS = 0x44
};

__attribute__ ((__nonnull__))
void bmm_defhead(struct bmm_head*);

// TODO Use functions as filters or something.
__attribute__ ((__nonnull__))
bool bmm_msg_get(struct bmm_head*, struct bmm_dem_state*);

__attribute__ ((__nonnull__))
void bmm_msg_put(struct bmm_head const*, struct bmm_dem_state const*);

#endif
