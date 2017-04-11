#ifndef BMM_MSG_H
/// Messaging protocol.
#define BMM_MSG_H

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

#include "ext.h"
#include "dem.h"
#include "size.h"

/// This enumeration specifies message size.
enum bmm_msg_width {
  BMM_MSG_WIDTH_NARROW,
  BMM_MSG_WIDTH_WIDE
};

/// This enumeration distinguishes
/// literal-terminated messages from size-prefixed messages.
enum bmm_msg_tag {
  BMM_MSG_TAG_LT,
  BMM_MSG_TAG_SP
};

/// This enumeration specifies integer endianness.
enum bmm_msg_endian {
  BMM_MSG_ENDIAN_LITTLE,
  BMM_MSG_ENDIAN_BIG
};

/// This structure allows the user to signal
/// what kind of message they want to send.
struct bmm_msg_choice {
  enum bmm_msg_width width;
  enum bmm_msg_endian endian;
  enum bmm_msg_tag tag;
  union {
    struct {
      size_t nmemb;
      uint8_t buf[8];
    } term;
    size_t size;
  } msg;
};

#define BMM_MSG_BIT_WIDE 7
#define BMM_MSG_BIT_VAR 3
#define BMM_MSG_BIT_TAG 2

/// The buffer `ptr` must be able to hold at least `*pnmemb` octets.
__attribute__ ((__nonnull__))
inline void bmm_msg_buf(uint8_t* const pbuf, size_t* const pnmemb,
    struct bmm_msg_choice const* const choice) {
  uint8_t userset0;

  switch (choice->endian) {
    case BMM_MSG_ENDIAN_LITTLE:
      userset0 = 0;

      break;
    case BMM_MSG_ENDIAN_BIG:
      userset0 = 7;

      break;
  }

  userset0 <<= 4;

  uint8_t derived0;

  uint8_t buf[8];
  size_t nmemb;

  switch (choice->tag) {
    case BMM_MSG_TAG_LT:
      {
        dynamic_assert(choice->msg.term.nmemb > sizeof choice->msg.term.buf,
            "Buffer would overflow");

        size_t const logsize = bmm_size_cilog(choice->msg.term.nmemb, 2);

        derived0 = (uint8_t) logsize;

        derived0 = BMM_SETBIT(derived0, BMM_MSG_BIT_VAR);

        (void) memcpy(buf, choice->msg.term.buf, choice->msg.term.nmemb);

        nmemb = choice->msg.term.nmemb;
      }

      break;
    case BMM_MSG_TAG_SP:
      {
        size_t const size = choice->msg.size;

        if (choice->msg.size < 8) {
          derived0 = (uint8_t) size;

          nmemb = 0;
        } else {
          size_t const logsize = bmm_size_cilog(size, 2);
          size_t const log2size = bmm_size_cilog(logsize, 2);
          size_t const powsize = bmm_size_pow(log2size, 2);

          derived0 = (uint8_t) log2size;

          derived0 = BMM_SETBIT(derived0, BMM_MSG_BIT_TAG);

          // TODO Yeah...
          switch (choice->endian) {
            case BMM_MSG_ENDIAN_LITTLE:
              for (size_t i = 0; i < powsize; ++i)
                buf[i] = (uint8_t) (size >> i * 8 & 0xff);

              break;
            case BMM_MSG_ENDIAN_BIG:
              for (size_t i = 0; i < powsize; ++i)
                buf[powsize - 1 - i] = (uint8_t) (size >> i * 8 & 0xff);

              break;
          }

          nmemb = powsize;
        }

        derived0 = BMM_SETBIT(derived0, BMM_MSG_BIT_VAR);
      }

      break;
  }

  size_t width0;

  size_t j = 1;

  // BMM_TESTBIT(userset0, BMM_MSG_BIT_WIDE)
  switch (choice->width) {
    case BMM_MSG_WIDTH_NARROW:
      width0 = 0;

      pbuf[0] = width0 | userset0 | derived0;

      ++j;

      break;
    case BMM_MSG_WIDTH_WIDE:
      width0 = 1 << 7;

      pbuf[0] = width0 | userset0;
      pbuf[1] = derived0;

      break;
  }

  for (size_t i = 0; i < nmemb; ++i)
    pbuf[j + i] = buf[i];

  *pnmemb = nmemb;
}

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
  BMM_MSG_NSTEP = 60,
  BMM_MSG_NPART = 142,
  BMM_MSG_PARTS = 144,
  BMM_MSG_NEIGH = 168,
  BMM_MSG_EKINE = 185,
  BMM_MSG_ESMOM = 186,
  BMM_MSG_EVMOM = 187
};

// TODO Wrangle prototypes.

__attribute__ ((__nonnull__ (2)))
bool bmm_msg_preread(size_t* const ptr, struct bmm_msg_head const* const head);

__attribute__ ((__nonnull__))
bool bmm_msg_prewrite(struct bmm_msg_head const* const bad_head,
    size_t const size);

__attribute__ ((__nonnull__))
void bmm_head_def(struct bmm_msg_head*);

__attribute__ ((__deprecated__, __nonnull__))
bool bmm_msg_get(struct bmm_msg_head*, struct bmm_dem*);

__attribute__ ((__deprecated__, __nonnull__))
void bmm_msg_put(struct bmm_msg_head const*, struct bmm_dem const*);

#endif
