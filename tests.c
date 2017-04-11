#include <cheat.h>
#include <cheats.h>
#include <stdint.h>
#include <string.h>

#include "msg.h"

CHEAT_DECLARE(
  static uint64_t const sample = 0xfadedacebeeffeed;

  struct bmm_msg_test {
    size_t i;
    uint8_t buf[10];
  };

  bool bmm_msg_test_read(uint8_t* const x, void* const ptr) {
    struct bmm_msg_test* const test = ptr;

    *x = test->buf[test->i];
    ++test->i;

    return true;
  }

  bool bmm_msg_test_write(uint8_t const* x, void* const ptr) {
    struct bmm_msg_test* const test = ptr;

    test->buf[test->i] = *x;
    ++test->i;

    return true;
  }
)

CHEAT_TEST(msg,
  for (enum bmm_msg_width width = BMM_MSG_WIDTH_NARROW;
      width < BMM_MSG_WIDTH_WIDE;
      ++width)
    for (enum bmm_msg_endian endian = BMM_MSG_ENDIAN_LITTLE;
        endian < BMM_MSG_ENDIAN_BIG;
        ++endian) {
      for (size_t e = 0; e < 4; ++e) {
        size_t const nmemb = 1 << e;

        struct bmm_msg_spec specin;
        struct bmm_msg_spec specout;

        specin.width = width;
        specin.endian = endian;
        specin.tag = BMM_MSG_TAG_LT;
        specin.msg.term.nmemb = nmemb;

        for (size_t i = 0; i < nmemb; ++i)
          specin.msg.term.buf[i] = sample << i * 8 & 0xff;

        struct bmm_msg_test test;

        test.i = 0;
        bmm_msg_spec_write(&specin, bmm_msg_test_write, &test);

        test.i = 0;
        bmm_msg_spec_read(&specout, bmm_msg_test_read, &test);

        cheat_assert_int(specin.width, specout.width);
        cheat_assert_int(specin.endian, specout.endian);
        cheat_assert_int(specin.tag, specout.tag);
        cheat_assert_size(specin.msg.term.nmemb, specout.msg.term.nmemb);

        for (size_t i = 0; i < specin.msg.term.nmemb; ++i)
          cheat_assert_unsigned_char(specin.msg.term.buf[i], specout.msg.term.buf[i]);
      }

      for (size_t size = 0; size < 300; ++size) {
        struct bmm_msg_spec specin;
        struct bmm_msg_spec specout;

        specin.width = width;
        specin.endian = endian;
        specin.tag = BMM_MSG_TAG_SP;
        specin.msg.size = size;

        struct bmm_msg_test test;

        test.i = 0;
        bmm_msg_spec_write(&specin, bmm_msg_test_write, &test);

        test.i = 0;
        bmm_msg_spec_read(&specout, bmm_msg_test_read, &test);

        cheat_assert_int(specin.width, specout.width);
        cheat_assert_int(specin.endian, specout.endian);
        cheat_assert_int(specin.tag, specout.tag);
        cheat_assert_size(specin.msg.size, specout.msg.size);
      }
    }
)
