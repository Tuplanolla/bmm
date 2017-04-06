#include <arpa/inet.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

static uint32_t const sample = 0xfadedace;

static uint32_t const lempsa = 0xcedadefa;

static bool bigend(void) {
  return htonl(sample) == sample;
}

static bool litend(void) {
  return htonl(sample) == lempsa;
}

int main(void) {
  return puts(bigend() ? "This system is in big endian." :
      litend() ? "This system is in little endian." :
      "This system is not in big or little endian.") == EOF ?
    EXIT_FAILURE : EXIT_SUCCESS;
}
