#include "clopts.h"
#include "errors.h"
#include "exts.h"
#include <stdbool.h>
#include <string.h>

bool bmm_clopts(char const* const* const args, size_t const narg,
    bool (* const f)(char const*, char const*, void*), void* const p) {
  char const* key = NULL;
  enum {KEY, VALUE} state = KEY;

  for (size_t iarg = 0; iarg < narg; ++iarg)
    switch (state) {
      case KEY:
        if (strncmp(args[iarg], "--", 2) != 0) {
          bmm_error("Invalid key '%s'.", args[iarg]);

          return false;
        } else {
          key = &args[iarg][2];
          state = VALUE;

          break;
        }
      case VALUE:
        if (!f(key, args[iarg], p)) {
          bmm_error("Invalid value '%s' for key '%s'.", args[iarg], key);

          return false;
        } else {
          state = KEY;

          break;
        }
    }

  switch (state) {
    case KEY:
      return true;
    case VALUE:
      bmm_error("No value for key '%s'.", args[narg - 1]);

      return false;
  }
}

