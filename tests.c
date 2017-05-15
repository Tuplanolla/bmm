#include <cheat.h>
#include <cheats.h>
#include <stdint.h>
#include <string.h>

#include "cpp.h"
#include "endy.h"
#include "ext.h"
// #include "moore.h"
#include "msg.h"
#include "size.h"

CHEAT_TEST(cpp_testbit,
  for (size_t i = 0; i < 8; ++i)
    cheat_assert(BMM_TESTBIT(1 << i, i));
)

CHEAT_TEST(size_hc,
  size_t const ndim = 3;
  size_t const nper = 3;

  size_t ij[3];

  for (size_t i = 0; i < 3 * 3 * 3; ++i) {
    bmm_size_hc(ij, i, ndim, nper);
    size_t const j = bmm_size_unhc(ij, ndim, nper);

    cheat_assert_size(j, i);
  }
)

CHEAT_TEST(size_hcd,
  size_t const ndim = 3;
  size_t const nper[] = {2, 3, 4};

  size_t ij[3];

  for (size_t i = 0; i < 2 * 3 * 4; ++i) {
    bmm_size_hcd(ij, i, ndim, nper);
    size_t const j = bmm_size_unhcd(ij, ndim, nper);

    cheat_assert_size(j, i);
  }
)

CHEAT_DECLARE(
  static enum bmm_msg_prio const msg_prio[] = {
    BMM_MSG_PRIO_LOW,
    BMM_MSG_PRIO_HIGH
  };

  static enum bmm_endy const msg_endy[] = {
    BMM_ENDY_LITTLE,
    BMM_ENDY_BIG
  };

  struct msg {
    size_t i;
    unsigned char buf[BMM_MSG_HEADSIZE];
  };

  static enum bmm_io_read msg_read(void* const pbuf, size_t const n,
      void* const ptr) {
    unsigned char* const buf = pbuf;
    struct msg* const msg = ptr;

    for (size_t i = 0; i < n; ++i) {
      buf[i] = msg->buf[msg->i];
      ++msg->i;
    }

    return BMM_IO_READ_SUCCESS;
  }

  static bool msg_write(void const* pbuf, size_t const n, void* const ptr) {
    unsigned char const* const buf = pbuf;
    struct msg* const msg = ptr;

    for (size_t i = 0; i < n; ++i) {
      msg->buf[msg->i] = buf[i];
      ++msg->i;
    }

    return true;
  }
)

CHEAT_TEST(msg_spec_sp_iso,
  for (size_t iprio = 0; iprio < nmembof(msg_prio); ++iprio)
    for (size_t iendy = 0; iendy < nmembof(msg_endy); ++iendy)
      for (size_t size = 0; size < 1000; ++size) {
        struct bmm_msg_spec out;

        out.prio = msg_prio[iprio];
        out.endy = msg_endy[iendy];
        out.tag = BMM_MSG_TAG_SP;
        out.msg.size = size;

        struct msg msg;

        msg.i = 0;
        cheat_assert(bmm_msg_spec_write(&out, msg_write, &msg));

        struct bmm_msg_spec in;

        msg.i = 0;
        cheat_assert(bmm_msg_spec_read(&in, msg_read, &msg) ==
            BMM_IO_READ_SUCCESS);

        cheat_assert_int(out.prio, in.prio);
        cheat_assert_int(out.endy, in.endy);
        cheat_assert_int(out.tag, in.tag);
        cheat_assert_size(out.msg.size, in.msg.size);
      }
)

CHEAT_TEST(msg_spec_lt_iso,
  for (size_t iprio = 0; iprio < nmembof(msg_prio); ++iprio)
    for (size_t iendy = 0; iendy < nmembof(msg_endy); ++iendy)
      for (size_t e = 0; e < 4; ++e) {
        struct bmm_msg_spec out;

        out.prio = msg_prio[iprio];
        out.endy = msg_endy[iendy];
        out.tag = BMM_MSG_TAG_LT;
        out.msg.term.e = e;

        for (size_t i = 0; i < bmm_size_pow(2, e); ++i)
          out.msg.term.buf[i] = i << i * 8 & 0xff;

        struct msg msg;

        msg.i = 0;
        cheat_assert(bmm_msg_spec_write(&out, msg_write, &msg));

        struct bmm_msg_spec in;

        msg.i = 0;
        cheat_assert(bmm_msg_spec_read(&in, msg_read, &msg) ==
            BMM_IO_READ_SUCCESS);

        cheat_assert_int(out.prio, in.prio);
        cheat_assert_int(out.endy, in.endy);
        cheat_assert_int(out.tag, in.tag);
        cheat_assert_size(out.msg.term.e, in.msg.term.e);

        for (size_t i = 0; i < out.msg.term.e; ++i)
          cheat_assert_unsigned_char(out.msg.term.buf[i], in.msg.term.buf[i]);
      }
)
