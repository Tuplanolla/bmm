#include <cheat.h>
#include <cheats.h>
#include <stdint.h>
#include <string.h>

#include "cpp.h"
#include "msg.h"
#include "size.h"

CHEAT_TEST(cpp_testbit,
  for (size_t i = 0; i < 8; ++i)
    cheat_assert(BMM_TESTBIT(1 << i, i));
)

CHEAT_DECLARE(
  static uint64_t const sample = 0xfadedacebeeffeed;

  struct msg_buf {
    size_t i;
    uint8_t buf[BMM_MSG_HEADSIZE];
  };

  static bool msg_read(uint8_t* const x, size_t const n, void* const ptr) {
    struct msg_buf* const buf = ptr;

    for (size_t i = 0; i < n; ++i) {
      x[i] = buf->buf[buf->i];
      ++buf->i;
    }

    return true;
  }

  static bool msg_write(uint8_t const* x, size_t const n, void* const ptr) {
    struct msg_buf* const buf = ptr;

    for (size_t i = 0; i < n; ++i) {
      buf->buf[buf->i] = x[i];
      ++buf->i;
    }

    return true;
  }
)

CHEAT_TEST(msg_spec_sp_iso,
  for (enum bmm_msg_prio prio = BMM_MSG_PRIO_LOW;
      prio < BMM_MSG_PRIO_HIGH;
      ++prio)
    for (enum bmm_msg_endian endian = BMM_MSG_ENDIAN_LITTLE;
        endian < BMM_MSG_ENDIAN_BIG;
        ++endian)
      for (size_t size = 0; size < 1000; ++size) {
        struct bmm_msg_spec in;

        in.prio = prio;
        in.endian = endian;
        in.tag = BMM_MSG_TAG_SP;
        in.msg.size = size;

        struct msg_buf buf;

        buf.i = 0;
        bmm_msg_spec_write(&in, msg_write, &buf);

        struct bmm_msg_spec out;

        buf.i = 0;
        bmm_msg_spec_read(&out, msg_read, &buf);

        cheat_assert_int(in.prio, out.prio);
        cheat_assert_int(in.endian, out.endian);
        cheat_assert_int(in.tag, out.tag);
        cheat_assert_size(in.msg.size, out.msg.size);
      }
)

CHEAT_TEST(msg_spec_lt_iso,
  for (enum bmm_msg_prio prio = BMM_MSG_PRIO_LOW;
      prio < BMM_MSG_PRIO_HIGH;
      ++prio)
    for (enum bmm_msg_endian endian = BMM_MSG_ENDIAN_LITTLE;
        endian < BMM_MSG_ENDIAN_BIG;
        ++endian)
      for (size_t e = 0; e < 4; ++e) {
        struct bmm_msg_spec in;

        in.prio = prio;
        in.endian = endian;
        in.tag = BMM_MSG_TAG_LT;
        in.msg.term.e = e;

        for (size_t i = 0; i < bmm_size_pow(2, e); ++i)
          in.msg.term.buf[i] = sample << i * 8 & 0xff;

        struct msg_buf buf;

        buf.i = 0;
        bmm_msg_spec_write(&in, msg_write, &buf);

        struct bmm_msg_spec out;

        buf.i = 0;
        bmm_msg_spec_read(&out, msg_read, &buf);

        cheat_assert_int(in.prio, out.prio);
        cheat_assert_int(in.endian, out.endian);
        cheat_assert_int(in.tag, out.tag);
        cheat_assert_size(in.msg.term.e, out.msg.term.e);

        for (size_t i = 0; i < in.msg.term.e; ++i)
          cheat_assert_unsigned_char(in.msg.term.buf[i], out.msg.term.buf[i]);
      }
)
