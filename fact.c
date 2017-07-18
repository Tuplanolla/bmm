#include <limits.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

static bool terms(size_t const n) {
  if (printf("(1") < 0)
    return false;

  for (size_t i = 0; i < n; ++i)
    if (printf(" * %zu", i + 1) < 0)
      return false;

  if (printf(")") < 0)
    return false;

  return true;
}

static bool dirs(size_t const n) {
  for (size_t i = 0; i <= n; ++i)
    if (printf("#define BMM_FACT_%zu() ", i) < 0 ||
        !terms(i) ||
        printf("\n") < 0)
      return false;

  return true;
}

int main(int const argc, char **const argv) {
  if (argc != 2)
    return EXIT_FAILURE;

  unsigned long int n = strtoul(argv[1], NULL, 10);
  if (n == ULONG_MAX)
    return EXIT_FAILURE;

  if (!dirs((size_t) n))
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
