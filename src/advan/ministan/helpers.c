#include "common.h"
#include <stdio.h>

/**
 * Prints an array of doubles with customizable formatting
 * @param arr: pointer to the array of doubles
 * @param size: number of elements in the array
 * @param precision: number of decimal places (default suggestion: 2)
 * @param separator: string to separate elements (e.g., ", ", " ", "\t")
 */
void print_double_array(const char *descr, const double *arr, size_t size) {
  printf("%s: ", descr);

  if (arr == NULL || size == 0) {
    printf("[empty]\n");
    return;
  }

  printf("[");
  for (size_t i = 0; i < size; i++) {
    printf("%.e", arr[i]);
    if (i < size - 1)
      printf(", ");
  }
  printf("]\n");
}

void print_double(const char *descr, const double number) {
  printf("%s: %e\n", descr, number);
}

void print_doubles(const char *descr, const double *arr, size_t size) {
  print_double_array(descr, arr, size);
}

void print_error(const char *descr) { printf("Error condition: %s\n", descr); }