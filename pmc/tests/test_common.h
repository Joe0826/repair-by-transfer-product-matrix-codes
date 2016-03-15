#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define talloc(type, num) (type *) malloc(sizeof(type)*(num))
#define mprint jerasure_print_matrix

// Print an array.
static void aprint(int* a, int len) {
  for(int i = 0; i < len; i++)
    printf("%d, ", a[i]);
  printf("\n");
}

// Print chunks (ie, node data).
static void cprint(char** chunk_ptrs, int nchunks, int chunk_size)
{
  for(int i = 0; i < nchunks; i++) {
    for(int j = 0; j < chunk_size; j++) {
      char c = chunk_ptrs[i][j];
      printf("%c", c);
    }
    printf("\n");
  }
  printf("\n");
}

// Messes with memory (allocs/deallocs),
// to help catch assumptions about uninitialized memory.
// Can be placed anywhere; should not affect behavior of VALID programs.
static void memory_monkey() {
    const int N = 17; // A prime.
    int* ptrs[N];
    int i;
    // Allocate.
    for (i = 0; i < N; i++) {
        size_t len = sizeof(int)*(200 + 10*i);
        ptrs[i] = (int*) malloc(len);
        memset(ptrs[i], 255, len);
    }
    // Free (out-of-order), re-allocate.
    for (int j = 0; i > 0; j = (j+3) % N, i--) {
        printf("i, j: %d, %d\n",  i, j);
        free(ptrs[j]);
        ptrs[j] = (int*) malloc(sizeof(int)*(130 + 13*i));
    }
    // Free all.
    for (i = 0; i < N; i++) {
        free(ptrs[i]);
    }
}
