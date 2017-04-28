#include <stdbool.h>
#include <stddef.h>
#include <string.h>

#include "opt.h"
#include "ext.h"
#include "tle.h"

bool bmm_opt_parse(char const* const* const args, size_t const narg,
    bool (* const f)(char const*, char const*, void*), void* const ptr) {
  char const* key = NULL;
  enum {KEY, VALUE} state = KEY;

  for (size_t iarg = 0; iarg < narg; ++iarg)
    switch (state) {
      case KEY:
        if (strncmp(args[iarg], "--", 2) != 0) {
          BMM_TLE_EXTS(BMM_TLE_NUM_PARSE, "Invalid key '%s'", args[iarg]);

          return false;
        } else {
          key = &args[iarg][2];
          state = VALUE;

          break;
        }
      case VALUE:
        if (!f(key, args[iarg], ptr)) {
          BMM_TLE_EXTS(BMM_TLE_NUM_PARSE, "Invalid key and value '%s' and '%s'",
              key, args[iarg]);

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
      BMM_TLE_EXTS(BMM_TLE_NUM_PARSE, "No value for key '%s'", args[narg - 1]);

      return false;
  }

  dynamic_assert(false, "Impossible");
}

