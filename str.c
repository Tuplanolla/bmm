#include <errno.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "str.h"
#include "tle.h"

bool bmm_str_strtou(unsigned int* const ptr, char const* const str) {
  char* endptr;
  errno = 0;
  long long int const x = strtoll(str, &endptr, 10);
  if (errno != 0 || endptr[0] != '\0' || x < 0) {
    BMM_TLE_EXTS(BMM_TLE_PARSE, "Cannot parse '%s' as an unsigned integer",
        str);

    return false;
  } else {
    if (ptr != NULL)
      *ptr = (unsigned int) x;

    return true;
  }
}

bool bmm_str_strtoz(size_t* const ptr, char const* const str) {
  char* endptr;
  errno = 0;
  long long int const x = strtoll(str, &endptr, 10);
  if (errno != 0 || endptr[0] != '\0' || x < 0) {
    BMM_TLE_EXTS(BMM_TLE_PARSE, "Cannot parse '%s' as a size", str);

    return false;
  } else {
    if (ptr != NULL)
      *ptr = (size_t) x;

    return true;
  }
}

bool bmm_str_strtod(double* const ptr, char const* const str) {
  char* endptr;
  errno = 0;
  double const x = strtod(str, &endptr);
  if (errno != 0 || endptr[0] != '\0') {
    BMM_TLE_EXTS(BMM_TLE_PARSE, "Cannot parse '%s' as a floating-point number",
        str);

    return false;
  } else {
    if (ptr != NULL)
      *ptr = (double) x;

    return true;
  }
}

bool bmm_str_strtob(bool* const ptr, char const* const str) {
  bool x;
  if (strcmp(str, "0") == 0 ||
      strcmp(str, "false") == 0 ||
      strcmp(str, "no") == 0)
    x = false;
  else if (strcmp(str, "1") == 0 ||
      strcmp(str, "true") == 0 ||
      strcmp(str, "yes") == 0)
    x = true;
  else {
    BMM_TLE_EXTS(BMM_TLE_PARSE, "Cannot parse '%s' as a truth value", str);

    return false;
  }

  if (ptr != NULL)
    *ptr = x;

  return true;
}
