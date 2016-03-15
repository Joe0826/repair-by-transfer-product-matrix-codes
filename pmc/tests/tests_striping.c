#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../jerasure-kmg/include/galois.h"
#include "reed_sol.h"
#include "jerasure.h"

#include "MBR_MSR.h"
#include "product_matrix.h"

#define talloc(type, num) (type *) malloc(sizeof(type)*(num))

#define mprint jerasure_print_matrix

void aprint(int* a, int len) {
  for(int i = 0; i < len; i++)
    printf("%d, ", a[i]);
  printf("\n");
}

void testSingleEncode()
{
  const int w = 16;

  const int n = 5;
  const int d = 4;
  const int a = 4;
  const int k = 3; //for M below:
  const int B = 9;

  int* M = make_m_MBR(k, d);
  int* psi = make_psi_MBR(n, k, d, w); // n x d
  int* G = generator_for_PMC(psi, M, n, d, a, B);
  make_systematic(G, n, a, B, w);


  char* msg = "hello, world1!!!!!hello, world2!!!!!hello, world3!!!!!hello, world4!!!!";
  int msg_size = strlen(msg) + 1;
  //printf("%d\n", msg_size);

  char** dst_ptrs = talloc(char*, a);
  int chunk_size = msg_size / B; //bytes stored in all stripes for a particular coded-symbol.
  for (int i = 0; i < a; i++)
    dst_ptrs[i] = talloc(char, chunk_size);


  encode_striped(msg, G, dst_ptrs, a, B, msg_size, w);


  int num_stripes = msg_size / (B * w/8);
  printf("Node 0 stores %d stripes of alpha=%d 2-byte-symbols each:\n", num_stripes, a);
  for(int i = 0; i < a; i++) {
    for(int j = 0; j < msg_size / B; j++) {
      char c = dst_ptrs[i][j];
      printf("%c", c);
    }
    printf("\n");
  }
  printf("\n");
}

void cprint(char** chunk_ptrs, int nchunks, int chunk_size)
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

int main(int argc, char **argv)
{

  const int w = 16;

  const int n = 5;
  const int d = 4;
  const int a = 4;
  const int k = 3; //for M below:
  const int B = 9;

  int* M = make_m_MBR(k, d);
  int* psi = make_psi_MBR(n, k, d, w); // n x d
  int* G = generator_for_PMC(psi, M, n, d, a, B);
  make_systematic(G, n, a, B, w);


  char* msg = "hello, world1!!!!!hello, world2!!!!!hello, world3!!!!!hello, world4!!!!";
  int msg_size = strlen(msg) + 1;
  //printf("%d\n", msg_size);

  char** all_dst_ptrs = talloc(char*, n*a);
  int chunk_size = msg_size / B; //bytes stored in all stripes for a particular coded-symbol.
  for (int i = 0; i < n*a; i++)
    all_dst_ptrs[i] = talloc(char, chunk_size);


  /* Encode all nodes. */
  char** dst_i_ptr = all_dst_ptrs;
  int* G_i = G;
  for (int i = 0; i < n; i++) {
    encode_striped(msg, G_i, dst_i_ptr, a, B, msg_size, w);
    dst_i_ptr += a;
    G_i += a*B;
  }


  int num_stripes = msg_size / (B * w/8);
  char** node_ptrs = all_dst_ptrs + a*0;
  printf("Node stores %d stripes of alpha=%d 2-byte-symbols each:\n", num_stripes, a);

  cprint(node_ptrs, a, chunk_size);
  printf("\n");


  /* Test DC. */
  int nodes[] = {0, 1, 2};
  int* D = make_data_collect_matrix(G, k, a, B, nodes, w);

  char** chunk_ptrs = all_dst_ptrs;
  char* dst = talloc(char, msg_size);
  data_collect_striped(chunk_ptrs, D, dst, k, a, B, msg_size, w);
  printf("collected message: \n");
  printf("%s", dst);
  printf("\n");


  /* Test Repair. */
  int f = 4;
  int* mu_f = make_mu_MBR(f, k, d, w);

  // Helper nodes send...
  char** helper_chunks = talloc(char*, d);
  char** data_i_ptr = all_dst_ptrs;
  for (int i = 0; i < d; i++) {
    // Helper-node i sends helper_chunk[i] to the failed node.
    helper_chunks[i] = talloc(char, chunk_size);
    helper_symbol_striped(data_i_ptr, mu_f, helper_chunks[i], a, B, msg_size, w);
    data_i_ptr += a;
  }

  //Failed node repairs...
  int helper_nodes[] = {0, 1, 2, 3};
  int* A = make_repair_matrix(G, d, a, B, f, mu_f, helper_nodes, w);

  char** repaired_data_ptrs = talloc(char*, a);
  for (int i = 0; i < a; i++) {
    repaired_data_ptrs[i] = talloc(char, chunk_size);
  }

  repair_striped(helper_chunks, A, repaired_data_ptrs, d, a, B, msg_size, w);

  printf("\nfailed node chunks:\n");
  cprint(all_dst_ptrs + a*f, a, chunk_size);

  printf("repaired chunks:\n");
  cprint(repaired_data_ptrs, a, chunk_size);
}
