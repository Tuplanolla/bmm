#ifndef BMM_ENDY_H
/// Integer endianness.
#define BMM_ENDY_H

#include <stdbool.h>
#include <stdint.h>

#include "ext.h"

/// This enumeration specifies integer endianness.
enum bmm_endy {
  BMM_ENDY_LITTLE,
  BMM_ENDY_MIDDLE,
  BMM_ENDY_BIG
};

/// The call `bmm_endy_get()`
/// returns the integer endianness of the system.
__attribute__ ((__const__))
inline enum bmm_endy bmm_endy_get(void) {
#ifdef __BYTE_ORDER__

#if defined __ORDER_LITTLE_ENDIAN__ && \
  __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
  return BMM_ENDY_LITTLE;
#elif defined __ORDER_BIG_ENDIAN__ && \
  __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
  return BMM_ENDY_BIG;
#else
  return BMM_ENDY_MIDDLE;
#endif

#else

  bool little = true;
  bool big = true;

  static uint64_t const word = 0xbadfeed4fadedace;
  uint8_t const *const buf = (uint8_t const *) &word;

  for (size_t i = 0; i < 8; ++i) {
    uint8_t const octet = word >> i * 8 & 0xff;

    if (buf[i] != octet)
      little = false;
    else if (buf[n - 1 - i] != octet)
      big = false;
  }

  if (little)
    return BMM_ENDY_LITTLE;
  else if (big)
    return BMM_ENDY_BIG;
  else
    return BMM_ENDY_MIDDLE;

#endif
}

#endif
