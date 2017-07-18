#include <limits.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

static bool pars(size_t const n) {
  for (size_t i = 0; i < n; ++i) {
    if (i != 0)
      if (printf(", ") < 0)
        return false;

    if (printf("x%zu", i) < 0)
      return false;
  }

  return true;
}

static bool terms(size_t const n) {
  if (printf("(0") < 0)
    return false;

  for (size_t i = 0; i < n; ++i)
    if (printf(" | (1 << (x%zu))", i) < 0)
      return false;

  if (printf(")") < 0)
    return false;

  return true;
}

static bool dirs(size_t const n) {
  for (size_t i = 0; i <= n; ++i)
    if (printf("#define BMM_MASKBITS_%zu(", i) < 0 ||
        !pars(i) ||
        printf(") ") < 0 ||
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
