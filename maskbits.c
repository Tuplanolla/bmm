#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

static int cheat_print_parameters(size_t const parameters) {
  size_t index;

  for (index = 0;
      index < parameters;
      ++index) {
    if (index != 0)
      (void )printf(", ");
    (void )printf("x%zu", index);
  }
  return 0;
}

static int cheat_print_parametersx(size_t const parameters) {
  size_t index;

  (void )printf("(0");
  for (index = 0;
      index < parameters;
      ++index) {
    (void )printf(" | (1 << (x%zu))", index);
  }
  (void )printf(")");
  return 0;
}

static int cheat_print_definitions(size_t const parameters) {
  size_t index;

  for (index = 0;
      index <= parameters;
      ++index) {
    (void )printf("#define BMM_MASKBITS_%zu(", index);
    (void )cheat_print_parameters(index);
    (void )printf(") ");
    (void )cheat_print_parametersx(index);
    (void )printf("\n");
  }
  return 0;
}

int main(int const count, char** const arguments) {
  long int index;

  if (count != 2)
    return EXIT_FAILURE;

  index = strtol(arguments[1], NULL, 10);
  if (index < 1 || index == LONG_MAX)
    return EXIT_FAILURE;

  (void )cheat_print_definitions((size_t )index);

  return EXIT_SUCCESS;
}
