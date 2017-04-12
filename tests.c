#include <cheat.h>
#include <cheats.h>
#include <stdint.h>
#include <string.h>

#include "msg.h"
#include "size.h"

CHEAT_DECLARE(
  static uint64_t const sample = 0xfadedacebeeffeed;

  struct msg {
    size_t i;
    uint8_t buf[BMM_MSG_HEADSIZE];
  };

  static bool msg_read(uint8_t* const buf, size_t const n, void* const ptr) {
    struct msg* const msg = ptr;

    for (size_t i = 0; i < n; ++i) {
      buf[i] = msg->buf[msg->i];
      ++msg->i;
    }

    return true;
  }

  static bool msg_write(uint8_t const* buf, size_t const n, void* const ptr) {
    struct msg* const msg = ptr;

    for (size_t i = 0; i < n; ++i) {
      msg->buf[msg->i] = buf[i];
      ++msg->i;
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
        struct bmm_msg_spec out;

        out.prio = prio;
        out.endian = endian;
        out.tag = BMM_MSG_TAG_SP;
        out.msg.size = size;

        struct msg buf;

        buf.i = 0;
        bmm_msg_spec_write(&out, msg_write, &buf);

        struct bmm_msg_spec in;

        buf.i = 0;
        bmm_msg_spec_read(&in, msg_read, &buf);

        cheat_assert_int(out.prio, in.prio);
        cheat_assert_int(out.endian, in.endian);
        cheat_assert_int(out.tag, in.tag);
        cheat_assert_size(out.msg.size, in.msg.size);
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
        struct bmm_msg_spec out;

        out.prio = prio;
        out.endian = endian;
        out.tag = BMM_MSG_TAG_LT;
        out.msg.term.e = e;

        for (size_t i = 0; i < bmm_size_pow(2, e); ++i)
          out.msg.term.buf[i] = sample << i * 8 & 0xff;

        struct msg msg;

        msg.i = 0;
        bmm_msg_spec_write(&out, msg_write, &msg);

        struct bmm_msg_spec in;

        msg.i = 0;
        bmm_msg_spec_read(&in, msg_read, &msg);

        cheat_assert_int(out.prio, in.prio);
        cheat_assert_int(out.endian, in.endian);
        cheat_assert_int(out.tag, in.tag);
        cheat_assert_size(out.msg.term.e, in.msg.term.e);

        for (size_t i = 0; i < out.msg.term.e; ++i)
          cheat_assert_unsigned_char(out.msg.term.buf[i], in.msg.term.buf[i]);
      }
)
