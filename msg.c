#include <limits.h>
#include <stdbool.h>
#include <stdint.h>

#include "cpp.h"
#include "endy.h"
#include "io.h"
#include "msg.h"
#include "size.h"
#include "tle.h"

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

enum bmm_io_read bmm_msg_type_read(enum bmm_msg_type* const type,
    bmm_msg_reader const f, void* const ptr) {
  unsigned char buf;
  switch (f(&buf, 1, ptr)) {
    case BMM_IO_READ_ERROR:
      return BMM_IO_READ_ERROR;
    case BMM_IO_READ_EOF:
      return BMM_IO_READ_EOF;
  }

  *type = (enum bmm_msg_type) buf;

  return BMM_IO_READ_SUCCESS;
}

bool bmm_msg_type_write(enum bmm_msg_type const* const type,
    bmm_msg_writer const f, void* const ptr) {
  dynamic_assert(*type <= 0xff, "Type would be truncated");

  unsigned char buf = (unsigned char) *type;
  if (!f(&buf, 1, ptr))
    return false;

  return true;
}
