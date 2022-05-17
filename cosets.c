#include <stdbool.h>
#include <stdio.h>

int main(int argc, char ** argv) {
  if(argc < 2) {
    printf("There are not enough arguments\n");
    return 1;
  }

  size_t N;
  sscanf(argv[1], "%lu", &N);
  size_t multgroup = (1L << N) - 1; /* size of multiplicative group */

  for(size_t i = 1; i < (1L << N) - 1; i += 2) {
    /* Generate coset and check if i is the smallest guy in there */
    _Bool can_be_smallest = true;
    size_t t = i;
    for(size_t k = 1; k < N; ++k) {
      t = (t * 2) % multgroup;
      if (t < i) {
	can_be_smallest = false;
	break;
      }
    }
    if(can_be_smallest) {
      printf("%lu, ", i);
    }
  }
  printf("\n");

  return 0;
}
