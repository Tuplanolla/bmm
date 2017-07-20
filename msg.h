#ifndef BMM_MSG_H
/// Messaging protocol.
#define BMM_MSG_H

#include <stdbool.h>

#include "cpp.h"
#include "endy.h"
#include "ext.h"
#include "io.h"

/// These preprocessor directives define
/// the maximal number of octets for various parts of a message.
#define BMM_MSG_FLAGSIZE 1
#define BMM_MSG_PRESIZE 8
#define BMM_MSG_HEADSIZE (BMM_MSG_FLAGSIZE + BMM_MSG_PRESIZE)
#define BMM_MSG_NUMSIZE 1

/// This enumeration specifies message priority.
enum bmm_msg_prio {
  BMM_MSG_PRIO_LOW,
  BMM_MSG_PRIO_HIGH
};

/// This enumeration distinguishes
/// size-prefixed messages from literal-terminated messages.
enum bmm_msg_tag {
  BMM_MSG_TAG_SP,
  BMM_MSG_TAG_LT
};

/// This structure allows the user to signal
/// what kind of message they want to send.
/// Middle-endianness or free patterns are not supported.
struct bmm_msg_spec {
  enum bmm_msg_prio prio;
  enum bmm_endy endy;
  union {
    size_t size;
    struct {
      size_t e;
      unsigned char buf[BMM_MSG_PRESIZE];
    } term;
  } msg;
  enum bmm_msg_tag tag;
};

/// These preprocessor directives help work with bits in message headers.
#define BMM_MSG_MASK_PRIO (BMM_MASKBITS(7))
#define BMM_MSG_MASK_ENDY (BMM_MASKBITS(6, 5, 4))
#define BMM_MSG_MASK_VAR (BMM_MASKBITS(3))
#define BMM_MSG_MASK_TAG (BMM_MASKBITS(2))
#define BMM_MSG_MASK_FIXSIZE (BMM_MASKBITS(2, 1, 0))
#define BMM_MSG_MASK_VARSIZE (BMM_MASKBITS(1, 0))

/// The call `bmm_msg_spec_def(spec)`
/// writes the default message specification into `spec`.
__attribute__ ((__nonnull__))
void bmm_msg_spec_def(struct bmm_msg_spec *);

/// Assuming `bmm_msg_reader f`, the call `f(buf, n, ptr)`
/// reads `n` bytes into the buffer `buf`.
/// The additional `ptr` can be used for passing in a closure.
typedef enum bmm_io_read (*bmm_msg_reader)(void *, size_t, void *);

/// Assuming `bmm_msg_writer f`, the call `f(buf, n, ptr)`
/// writes `n` bytes from the buffer `buf`.
/// The additional `ptr` can be used for passing in a closure.
typedef bool (*bmm_msg_writer)(void const *, size_t, void *);

/// The call `bmm_msg_spec_read(spec, f, ptr)`
/// extracts the message specification `spec`
/// from the message header `buf` of length `n`
/// that is read by calling `f(buf, n, ptr)` as shown or
/// in chunks of size `0 < k < n`.
/// It is guaranteed that `n <= BMM_MSG_HEADSIZE`.
__attribute__ ((__nonnull__ (1, 2)))
enum bmm_io_read bmm_msg_spec_read(struct bmm_msg_spec *,
    bmm_msg_reader, void *);

/// The call `bmm_msg_spec_write(spec, f, ptr)`
/// builds the message header `buf` of length `n`
/// for the message specification `spec` and
/// writes it by sequentially calling `f(buf, n, ptr)` as shown or
/// in chunks of size `0 < k < n`.
/// It is guaranteed that `n <= BMM_MSG_HEADSIZE`.
__attribute__ ((__nonnull__ (1, 2)))
bool bmm_msg_spec_write(struct bmm_msg_spec const *, bmm_msg_writer, void *);

/// This enumeration is generated for message numbers.
enum bmm_msg_num {
#define BMM_MSG_DECLARE(id, num) \
  BMM_MSG_NUM_##id = num,
#include "msg_.h"
#undef BMM_MSG_DECLARE
};

/// The call `bmm_msg_to_str(ptr, num)`
/// sets the pointer `ptr` to
/// the string representation of the message number `num`.
bool bmm_msg_to_str(char const **, enum bmm_msg_num);

/// The call `bmm_msg_from_str(ptr, str)`
/// sets the pointer `ptr` to
/// the message number of the string representation `str`.
bool bmm_msg_from_str(enum bmm_msg_num *, char const *);

/// The call `bmm_msg_num_read(num, f, ptr)`
/// extracts the message number `num`
/// from the message header `buf` of length `n`
/// that is read by calling `f(buf, n, ptr)` as shown or
/// in chunks of size `0 < k < n`.
/// It is guaranteed that `n <= BMM_MSG_NUMSIZE`.
__attribute__ ((__nonnull__ (1, 2)))
enum bmm_io_read bmm_msg_num_read(enum bmm_msg_num *, bmm_msg_reader, void *);

/// The call `bmm_msg_num_write(num, f, ptr)`
/// builds the message header `buf` of length `n`
/// for the message number `num` and
/// writes it by sequentially calling `f(buf, n, ptr)` as shown or
/// in chunks of size `0 < k < n`.
/// It is guaranteed that `n <= BMM_MSG_NUMSIZE`.
__attribute__ ((__nonnull__ (1, 2)))
bool bmm_msg_num_write(enum bmm_msg_num const *, bmm_msg_writer, void *);

#endif
