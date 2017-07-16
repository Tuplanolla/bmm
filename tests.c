#include <cheat.h>
#include <cheats.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

#include "cpp.h"
#include "endy.h"
#include "ext.h"
#include "fp.h"
#include "geom2d.h"
#include "neigh.h"
#include "msg.h"
#include "size.h"

CHEAT_DECLARE(
  static size_t const ndim = 2;
  static size_t const nper[] = {6, 5};
  static bool const per[] = {true, false};
)

CHEAT_TEST(size_fact,
  cheat_assert_size(bmm_size_fact(0, 1), 1);
  cheat_assert_size(bmm_size_fact(1, 1), 1);
  cheat_assert_size(bmm_size_fact(2, 1), 2);
  cheat_assert_size(bmm_size_fact(3, 1), 6);
  cheat_assert_size(bmm_size_fact(4, 1), 24);
  cheat_assert_size(bmm_size_fact(5, 1), 120);
)

CHEAT_TEST(size_fact2,
  cheat_assert_size(bmm_size_fact(0, 2), 1);
  cheat_assert_size(bmm_size_fact(1, 2), 1);
  cheat_assert_size(bmm_size_fact(2, 2), 2);
  cheat_assert_size(bmm_size_fact(3, 2), 3);
  cheat_assert_size(bmm_size_fact(4, 2), 8);
  cheat_assert_size(bmm_size_fact(5, 2), 15);
)

CHEAT_TEST(size_flog,
  cheat_assert_size(bmm_size_flog(1, 2), 0);
  cheat_assert_size(bmm_size_flog(2, 2), 1);
  cheat_assert_size(bmm_size_flog(3, 2), 1);
  cheat_assert_size(bmm_size_flog(4, 2), 2);
  cheat_assert_size(bmm_size_flog(5, 2), 2);
  cheat_assert_size(bmm_size_flog(6, 2), 2);
)

CHEAT_TEST(size_clog,
  cheat_assert_size(bmm_size_clog(1, 2), 0);
  cheat_assert_size(bmm_size_clog(2, 2), 1);
  cheat_assert_size(bmm_size_clog(3, 2), 2);
  cheat_assert_size(bmm_size_clog(4, 2), 2);
  cheat_assert_size(bmm_size_clog(5, 2), 3);
  cheat_assert_size(bmm_size_clog(6, 2), 3);
)

CHEAT_TEST(geom2d_shell_inside,
  double const x[] = {0.5, 0.5};
  double const r = 1.0 / sqrt(3.0);
  double const xper[] = {1.0, 1.0};
  bool const per[] = {false, false};

  cheat_assert_double(bmm_geom2d_shellvol(x, r, xper, per),
      (M_2PI / 3.0) * r, 1e-6);
)

CHEAT_TEST(geom2d_shell_outside,
  double const x[] = {-0.5, -0.5};
  double const r = 1.0;
  double const xper[] = {1.0, 1.0};
  bool const per[] = {false, false};

  cheat_assert_double(bmm_geom2d_shellvol(x, r, xper, per),
      (M_PI_2 / 3.0) * r, 1e-6);
)

CHEAT_TEST(size_hc_ord,
  size_t ij[2];

  bmm_size_hc(ij, 0, ndim, nper[1]);
  cheat_assert_size(ij[0], 0);
  cheat_assert_size(ij[1], 0);

  bmm_size_hc(ij, 1, ndim, nper[1]);
  cheat_assert_size(ij[0], 0);
  cheat_assert_size(ij[1], 1);

  bmm_size_hc(ij, 2, ndim, nper[1]);
  cheat_assert_size(ij[0], 0);
  cheat_assert_size(ij[1], 2);
)

CHEAT_TEST(size_hc_iso,
  size_t ij[2];

  for (size_t i = 0; i < bmm_size_pow(nper[1], ndim); ++i) {
    bmm_size_hc(ij, i, ndim, nper[1]);
    size_t const j = bmm_size_unhc(ij, ndim, nper[1]);

    cheat_assert_size(j, i);
  }
)

CHEAT_TEST(size_hcd_ord,
  size_t ij[2];

  bmm_size_hcd(ij, 0, ndim, nper);
  cheat_assert_size(ij[0], 0);
  cheat_assert_size(ij[1], 0);

  bmm_size_hcd(ij, 1, ndim, nper);
  cheat_assert_size(ij[0], 0);
  cheat_assert_size(ij[1], 1);

  bmm_size_hcd(ij, 2, ndim, nper);
  cheat_assert_size(ij[0], 0);
  cheat_assert_size(ij[1], 2);
)

CHEAT_TEST(size_hcd_iso,
  size_t ij[2];

  for (size_t i = 0; i < bmm_size_prod(nper, ndim); ++i) {
    bmm_size_hcd(ij, i, ndim, nper);
    size_t const j = bmm_size_unhcd(ij, ndim, nper);

    cheat_assert_size(j, i);
  }
)

/*
CHEAT_TEST(neigh_np,
  cheat_assert_size(bmm_neigh_np(ndim), 9);
)

CHEAT_TEST(neigh_n,
  size_t buf[2];

  cheat_assert_size(bmm_neigh_n(buf, (size_t const[]) {0, 0}, ndim, nper), 4);
  cheat_assert_size(bmm_neigh_n(buf, (size_t const[]) {0, 1}, ndim, nper), 6);
  cheat_assert_size(bmm_neigh_n(buf, (size_t const[]) {1, 0}, ndim, nper), 6);
  cheat_assert_size(bmm_neigh_n(buf, (size_t const[]) {1, 1}, ndim, nper), 9);
)
*/

CHEAT_TEST(neigh_ncp,
  int const mask = BMM_NEIGH_MASK_FULL;

  cheat_assert_size(bmm_neigh_ncpij((size_t const[]) {1, 0},
        ndim, nper, per, mask), 6);
  cheat_assert_size(bmm_neigh_ncpij((size_t const[]) {1, 1},
        ndim, nper, per, mask), 9);
  cheat_assert_size(bmm_neigh_ncpij((size_t const[]) {1, 2},
        ndim, nper, per, mask), 9);
  cheat_assert_size(bmm_neigh_ncpij((size_t const[]) {2, 0},
        ndim, nper, per, mask), 6);
  cheat_assert_size(bmm_neigh_ncpij((size_t const[]) {2, 1},
        ndim, nper, per, mask), 9);
  cheat_assert_size(bmm_neigh_ncpij((size_t const[]) {2, 2},
        ndim, nper, per, mask), 9);
  cheat_assert_size(bmm_neigh_ncpij((size_t const[]) {3, 0},
        ndim, nper, per, mask), 6);
  cheat_assert_size(bmm_neigh_ncpij((size_t const[]) {3, 1},
        ndim, nper, per, mask), 9);
  cheat_assert_size(bmm_neigh_ncpij((size_t const[]) {3, 2},
        ndim, nper, per, mask), 9);
)

/*
CHEAT_TEST(neigh_ijp,
  size_t ij[2];

  bmm_neigh_ijp(ij, (size_t const[]) {1, 1}, 0, ndim, nper);
  cheat_assert_size(ij[0], 4);
  cheat_assert_size(ij[1], 3);

  bmm_neigh_ijp(ij, (size_t const[]) {1, 1}, 1, ndim, nper);
  cheat_assert_size(ij[0], 4);
  cheat_assert_size(ij[1], 1);

  bmm_neigh_ijp(ij, (size_t const[]) {1, 1}, 2, ndim, nper);
  cheat_assert_size(ij[0], 4);
  cheat_assert_size(ij[1], 2);

  bmm_neigh_ijp(ij, (size_t const[]) {1, 1}, 3, ndim, nper);
  cheat_assert_size(ij[0], 1);
  cheat_assert_size(ij[1], 3);

  bmm_neigh_ijp(ij, (size_t const[]) {1, 1}, 4, ndim, nper);
  cheat_assert_size(ij[0], 1);
  cheat_assert_size(ij[1], 1);

  bmm_neigh_ijp(ij, (size_t const[]) {1, 1}, 5, ndim, nper);
  cheat_assert_size(ij[0], 1);
  cheat_assert_size(ij[1], 2);

  bmm_neigh_ijp(ij, (size_t const[]) {1, 1}, 6, ndim, nper);
  cheat_assert_size(ij[0], 2);
  cheat_assert_size(ij[1], 3);

  bmm_neigh_ijp(ij, (size_t const[]) {1, 1}, 7, ndim, nper);
  cheat_assert_size(ij[0], 2);
  cheat_assert_size(ij[1], 1);

  bmm_neigh_ijp(ij, (size_t const[]) {1, 1}, 8, ndim, nper);
  cheat_assert_size(ij[0], 2);
  cheat_assert_size(ij[1], 2);
)

CHEAT_TEST(neigh_ij,
  size_t ij[2];

  bmm_neigh_ij(ij, (size_t const[]) {1, 1}, 0, ndim, nper);
  cheat_assert_size(ij[0], 0);
  cheat_assert_size(ij[1], 0);

  bmm_neigh_ij(ij, (size_t const[]) {1, 1}, 1, ndim, nper);
  cheat_assert_size(ij[0], 0);
  cheat_assert_size(ij[1], 1);

  bmm_neigh_ij(ij, (size_t const[]) {1, 1}, 2, ndim, nper);
  cheat_assert_size(ij[0], 0);
  cheat_assert_size(ij[1], 2);

  bmm_neigh_ij(ij, (size_t const[]) {1, 1}, 3, ndim, nper);
  cheat_assert_size(ij[0], 1);
  cheat_assert_size(ij[1], 0);

  bmm_neigh_ij(ij, (size_t const[]) {1, 1}, 4, ndim, nper);
  cheat_assert_size(ij[0], 1);
  cheat_assert_size(ij[1], 1);

  bmm_neigh_ij(ij, (size_t const[]) {1, 1}, 5, ndim, nper);
  cheat_assert_size(ij[0], 1);
  cheat_assert_size(ij[1], 2);

  bmm_neigh_ij(ij, (size_t const[]) {1, 1}, 6, ndim, nper);
  cheat_assert_size(ij[0], 2);
  cheat_assert_size(ij[1], 0);

  bmm_neigh_ij(ij, (size_t const[]) {1, 1}, 7, ndim, nper);
  cheat_assert_size(ij[0], 2);
  cheat_assert_size(ij[1], 1);

  bmm_neigh_ij(ij, (size_t const[]) {1, 1}, 8, ndim, nper);
  cheat_assert_size(ij[0], 2);
  cheat_assert_size(ij[1], 2);
)
*/

CHEAT_TEST(neigh_ijcp,
  int const mask = BMM_NEIGH_MASK_FULL;

  size_t ij[2];

  bmm_neigh_ijcpij(ij, (size_t const[]) {1, 1}, 0, ndim, nper, per, mask);
  cheat_assert_size(ij[0], 4);
  cheat_assert_size(ij[1], 0);

  bmm_neigh_ijcpij(ij, (size_t const[]) {1, 1}, 1, ndim, nper, per, mask);
  cheat_assert_size(ij[0], 4);
  cheat_assert_size(ij[1], 1);

  bmm_neigh_ijcpij(ij, (size_t const[]) {1, 1}, 2, ndim, nper, per, mask);
  cheat_assert_size(ij[0], 4);
  cheat_assert_size(ij[1], 2);

  bmm_neigh_ijcpij(ij, (size_t const[]) {1, 1}, 3, ndim, nper, per, mask);
  cheat_assert_size(ij[0], 1);
  cheat_assert_size(ij[1], 0);

  bmm_neigh_ijcpij(ij, (size_t const[]) {1, 1}, 4, ndim, nper, per, mask);
  cheat_assert_size(ij[0], 1);
  cheat_assert_size(ij[1], 1);

  bmm_neigh_ijcpij(ij, (size_t const[]) {1, 1}, 5, ndim, nper, per, mask);
  cheat_assert_size(ij[0], 1);
  cheat_assert_size(ij[1], 2);

  bmm_neigh_ijcpij(ij, (size_t const[]) {1, 1}, 6, ndim, nper, per, mask);
  cheat_assert_size(ij[0], 2);
  cheat_assert_size(ij[1], 0);

  bmm_neigh_ijcpij(ij, (size_t const[]) {1, 1}, 7, ndim, nper, per, mask);
  cheat_assert_size(ij[0], 2);
  cheat_assert_size(ij[1], 1);

  bmm_neigh_ijcpij(ij, (size_t const[]) {1, 1}, 8, ndim, nper, per, mask);
  cheat_assert_size(ij[0], 2);
  cheat_assert_size(ij[1], 2);
)

/*
CHEAT_TEST(neigh_ijpuh,
  size_t ij[2];

  bmm_neigh_ijpuh(ij, (size_t const[]) {4, 4}, 0, ndim, nper);
  cheat_assert_size(ij[0], 4);
  cheat_assert_size(ij[1], 1);

  bmm_neigh_ijpuh(ij, (size_t const[]) {4, 4}, 1, ndim, nper);
  cheat_assert_size(ij[0], 4);
  cheat_assert_size(ij[1], 2);

  bmm_neigh_ijpuh(ij, (size_t const[]) {4, 4}, 2, ndim, nper);
  cheat_assert_size(ij[0], 1);
  cheat_assert_size(ij[1], 3);
)

CHEAT_TEST(neigh_ijuh,
  size_t ij[2];

  bmm_neigh_ijuh(ij, (size_t const[]) {4, 4}, 0, ndim, nper);
  cheat_assert_size(ij[0], 4);
  cheat_assert_size(ij[1], 4);

  bmm_neigh_ijuh(ij, (size_t const[]) {4, 4}, 1, ndim, nper);
  cheat_assert_size(ij[0], 5);
  cheat_assert_size(ij[1], 3);

  bmm_neigh_ijuh(ij, (size_t const[]) {4, 4}, 2, ndim, nper);
  cheat_assert_size(ij[0], 5);
  cheat_assert_size(ij[1], 4);
)
*/

CHEAT_TEST(neigh_ijcpuh,
  int const mask = BMM_NEIGH_MASK_UPPERH;

  size_t ij[2];

  bmm_neigh_ijcpij(ij, (size_t const[]) {4, 4}, 0, ndim, nper, per, mask);
  cheat_assert_size(ij[0], 4);
  cheat_assert_size(ij[1], 4);

  bmm_neigh_ijcpij(ij, (size_t const[]) {4, 4}, 1, ndim, nper, per, mask);
  cheat_assert_size(ij[0], 1);
  cheat_assert_size(ij[1], 3);

  bmm_neigh_ijcpij(ij, (size_t const[]) {4, 4}, 2, ndim, nper, per, mask);
  cheat_assert_size(ij[0], 1);
  cheat_assert_size(ij[1], 4);
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
