#include <errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

static void moi(long double* const meanptr,
    long double* const varptr,
    long double* const sdptr,
    long double* const semptr,
    long long unsigned int const dims,
    long long unsigned int const iters) {
  long double m1 = 0.0;
  long double m2 = 0.0;

  long long unsigned int inball = 0;

  gsl_rng* const rng = gsl_rng_alloc(gsl_rng_mt19937);

  for (long long unsigned int incube = 0; incube < iters; ++incube) {
    long double proj2 = 0.0;

    for (long long unsigned int axis = 1; axis < dims; ++axis)
      proj2 += (long double) gsl_pow_2(gsl_rng_uniform(rng));

    long double const diag2 = proj2 +
      (long double) gsl_pow_2(gsl_rng_uniform(rng));

    if (diag2 < 1.0) {
      ++inball;

      long double const dproj2 = proj2 - m1;

      m1 += dproj2 / (long double) inball;
      m2 += dproj2 * (proj2 - m1);
    }
  }

  gsl_rng_free(rng);

  long double const mean = inball == 0 ? NAN : m1;
  long double const var = inball == 0 ? NAN : inball == 1 ? 0.0 :
    m2 / (long double) (inball - 1);
  long double const sd = inball == 0 ? NAN : inball == 1 ? 0.0 : sqrt(var);
  long double const sem = inball == 0 ? NAN : inball == 1 ? 0.0 :
    sqrt(var / (long double) inball);

  if (meanptr != NULL)
    *meanptr = mean;
  if (varptr != NULL)
    *varptr = var;
  if (sdptr != NULL)
    *sdptr = sd;
  if (semptr != NULL)
    *semptr = sem;
}

static bool parse(long long unsigned int* const ptr, char const* const str) {
  char* endptr;
  errno = 0;
  long long unsigned int x = strtoull(str, &endptr, 10);
  if (errno != 0)
    return false;
  if (endptr[0] != '\0') {
    errno = EINVAL;

    return false;
  }

  if (ptr != NULL)
    *ptr = (long long unsigned int) x;

  return true;
}

int main(int const argc, char** const argv) {
  if (argc != 3) {
    errno = EINVAL;
    perror("main");

    return EXIT_FAILURE;
  }

  long long unsigned int dims;
  long long unsigned int iters;

  if (!parse(&dims, argv[1]) || !parse(&iters, argv[2])) {
    perror("parse");

    return EXIT_FAILURE;
  }

  long double mean;
  long double sem;
  moi(&mean, NULL, NULL, &sem, dims, iters);

  if (printf("d = %llu, n = %llu, j = %Lg +- %Lg\n",
        dims, iters, mean, sem) < 0) {
    perror("printf");

    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
