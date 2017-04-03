#include "opt.h"
#include "err.h"
#include <stdbool.h>
#include <stddef.h>
#include <string.h>

bool bmm_opt_parse(char const* const* const args, size_t const narg,
    bool (* const f)(char const*, char const*, void*), void* const ptr) {
  char const* key = NULL;
  enum {KEY, VALUE} state = KEY;

  for (size_t iarg = 0; iarg < narg; ++iarg)
    switch (state) {
      case KEY:
        if (strncmp(args[iarg], "--", 2) != 0) {
          BMM_ERR_FWARN(NULL, "Invalid key '%s'", args[iarg]);

          return false;
        } else {
          key = &args[iarg][2];
          state = VALUE;

          break;
        }
      case VALUE:
        if (!f(key, args[iarg], ptr)) {
          BMM_ERR_FWARN(NULL, "Invalid key and value '%s' and '%s'",
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
      BMM_ERR_FWARN(NULL, "No value for key '%s'", args[narg - 1]);

      return false;
  }

  BMM_ERR_ABORT(NULL);
}

