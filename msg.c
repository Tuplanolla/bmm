#include <ctype.h>
#include <limits.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "conf.h"
#include "cpp.h"
#include "endy.h"
#include "ext.h"
#include "io.h"
#include "msg.h"
#include "size.h"
#include "tle.h"

static struct {
  enum bmm_msg_num num;
  char str[BMM_NCHARID];
  bool lowered;
} bmm_msg_decl[] = {
#define BMM_MSG_DECLARE(id, _) \
  {.num = BMM_MSG_NUM_##id, .str = #id, .lowered = false},
#include "msgnum.h"
#undef BMM_MSG_DECLARE
};

__attribute__ ((__pure__))
static size_t bmm_msg_find(enum bmm_msg_num const num) {
  for (size_t i = 0; i < sizeof bmm_msg_decl / sizeof *bmm_msg_decl; ++i)
    if (bmm_msg_decl[i].num == num)
      return i;

  return SIZE_MAX;
}

static void bmm_msg_lower(size_t const i) {
  if (!bmm_msg_decl[i].lowered) {
    for (size_t j = 0; bmm_msg_decl[i].str[j] != '\0'; ++j)
      bmm_msg_decl[i].str[j] = (char) tolower(bmm_msg_decl[i].str[j]);

    bmm_msg_decl[i].lowered = true;
  }
}

bool bmm_msg_to_str(char const** const ptr, enum bmm_msg_num const num) {
  size_t const i = bmm_msg_find(num);
  if (i == SIZE_MAX)
    return false;

  bmm_msg_lower(i);

  if (ptr != NULL)
    *ptr = bmm_msg_decl[i].str;

  return true;
}

bool bmm_msg_from_str(enum bmm_msg_num* const ptr, char const* const str) {
  for (size_t i = 0; i < sizeof bmm_msg_decl / sizeof *bmm_msg_decl; ++i) {
    bmm_msg_lower(i);

    if (strcmp(bmm_msg_decl[i].str, str) == 0) {
      if (ptr != NULL)
        *ptr = bmm_msg_decl[i].num;

      return true;
    }
  }

  return false;
}

void bmm_msg_spec_def(struct bmm_msg_spec* const spec) {
  spec->prio = BMM_MSG_PRIO_LOW;
  spec->endy = bmm_endy_get();
  spec->tag = BMM_MSG_TAG_SP;
  spec->msg.size = 0;
}

enum bmm_io_read bmm_msg_spec_read(struct bmm_msg_spec* const spec,
    bmm_msg_reader const f, void* const ptr) {
  unsigned char flags;
  switch (f(&flags, 1, ptr)) {
    case BMM_IO_READ_ERROR:
      return BMM_IO_READ_ERROR;
    case BMM_IO_READ_EOF:
      return BMM_IO_READ_EOF;
  }

  if (BMM_MASKALL(flags, BMM_MSG_MASK_PRIO))
    spec->prio = BMM_MSG_PRIO_HIGH;
  else
    spec->prio = BMM_MSG_PRIO_LOW;

  switch (flags & BMM_MSG_MASK_ENDY) {
    case BMM_MASKBITS_0():
      spec->endy = BMM_ENDY_LITTLE;

      break;
    case ~BMM_MASKBITS_0() & BMM_MSG_MASK_ENDY:
      spec->endy = BMM_ENDY_BIG;

      break;
    default:
      BMM_TLE_EXTS(BMM_TLE_UNIMPL, "Unsupported endianness");

      return BMM_IO_READ_ERROR;
  }

  if (BMM_MASKALL(flags, BMM_MSG_MASK_VAR)) {
    size_t const flagsize = (size_t) (flags & BMM_MSG_MASK_VARSIZE);
    size_t const presize = bmm_size_pow(2, flagsize);

    unsigned char buf[BMM_MSG_PRESIZE];
    switch (f(buf, presize, ptr)) {
      case BMM_IO_READ_ERROR:
        return BMM_IO_READ_ERROR;
      case BMM_IO_READ_EOF:
        return BMM_IO_READ_EOF;
    }

    if (BMM_MASKALL(flags, BMM_MSG_MASK_TAG)) {
      spec->tag = BMM_MSG_TAG_LT;

      spec->msg.term.e = flagsize;

      for (size_t i = 0; i < presize; ++i)
        spec->msg.term.buf[i] = buf[i];
    } else {
      spec->tag = BMM_MSG_TAG_SP;

      size_t size = 0;
      switch (spec->endy) {
        case BMM_ENDY_LITTLE:
          for (size_t i = 0; i < presize; ++i)
            size |= (size_t) buf[i] << i * CHAR_BIT;

          break;
        case BMM_ENDY_BIG:
          for (size_t i = 0; i < presize; ++i)
            size |= (size_t) buf[presize - 1 - i] << i * CHAR_BIT;

          break;
      }

      spec->msg.size = size;
    }
  } else {
    spec->tag = BMM_MSG_TAG_SP;

    spec->msg.size = (size_t) (flags & BMM_MSG_MASK_FIXSIZE);
  }

  return BMM_IO_READ_SUCCESS;
}

bool bmm_msg_spec_write(struct bmm_msg_spec const* const spec,
    bmm_msg_writer const f, void* const ptr) {
  unsigned char flags = BMM_MASKBITS_0();

  switch (spec->prio) {
    case BMM_MSG_PRIO_LOW:
      break;
    case BMM_MSG_PRIO_HIGH:
      flags |= BMM_MSG_MASK_PRIO;

      break;
  }

  switch (spec->endy) {
    case BMM_ENDY_LITTLE:
      break;
    case BMM_ENDY_MIDDLE:
      BMM_TLE_EXTS(BMM_TLE_UNIMPL, "Unsupported endianness");

      return false;
    case BMM_ENDY_BIG:
      flags |= ~BMM_MASKBITS_0() & BMM_MSG_MASK_ENDY;

      break;
  }

  size_t presize;
  unsigned char buf[BMM_MSG_PRESIZE];
  switch (spec->tag) {
    case BMM_MSG_TAG_SP:
      {
        size_t const size = spec->msg.size;

        if (size < bmm_size_pow(2, BMM_MSG_BITS_FIXSIZE)) {
          flags |= (unsigned char) size & BMM_MSG_MASK_FIXSIZE;

          presize = 0;
        } else {
          size_t const flagsize =
            bmm_size_clog(bmm_size_clog(size, CHAR_BIT), 2);

          flags |= BMM_MSG_MASK_VAR;
          flags |= (unsigned char) flagsize & BMM_MSG_MASK_VARSIZE;

          presize = bmm_size_pow(2, flagsize);

          switch (spec->endy) {
            case BMM_ENDY_LITTLE:
              for (size_t i = 0; i < presize; ++i)
                buf[i] = (unsigned char) (size >> i * CHAR_BIT & 0xff);

              break;
            case BMM_ENDY_BIG:
              for (size_t i = 0; i < presize; ++i)
                buf[presize - 1 - i] = (unsigned char)
                  (size >> i * CHAR_BIT & 0xff);

              break;
          }
        }
      }

      break;
    case BMM_MSG_TAG_LT:
      {
        size_t const flagsize = spec->msg.term.e;

        flags |= BMM_MSG_MASK_TAG;
        flags |= BMM_MSG_MASK_VAR;
        flags |= (unsigned char) flagsize & BMM_MSG_MASK_VARSIZE;

        presize = bmm_size_pow(2, flagsize);

        dynamic_assert(presize <= BMM_MSG_PRESIZE, "Buffer would overflow");

        for (size_t i = 0; i < presize; ++i)
          buf[i] = spec->msg.term.buf[i];
      }

      break;
  }

  if (!f(&flags, 1, ptr))
    return false;

  if (!f(buf, presize, ptr))
    return false;

  return true;
}

enum bmm_io_read bmm_msg_num_read(enum bmm_msg_num* const num,
    bmm_msg_reader const f, void* const ptr) {
  unsigned char buf;
  switch (f(&buf, 1, ptr)) {
    case BMM_IO_READ_ERROR:
      return BMM_IO_READ_ERROR;
    case BMM_IO_READ_EOF:
      return BMM_IO_READ_EOF;
  }

  *num = (enum bmm_msg_num) buf;

  return BMM_IO_READ_SUCCESS;
}

bool bmm_msg_num_write(enum bmm_msg_num const* const num,
    bmm_msg_writer const f, void* const ptr) {
  dynamic_assert(*num <= 0xff, "Type would be truncated");

  unsigned char buf = (unsigned char) *num;
  if (!f(&buf, 1, ptr))
    return false;

  return true;
}
