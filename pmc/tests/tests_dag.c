#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../jerasure-kmg/include/jerasure.h"

#include "product_matrix.h"
#include "dag_helper.h"
#include "MBR_MSR.h"

#define talloc(type, num) (type *) malloc(sizeof(type)*(num))

#define mprint jerasure_print_matrix

void aprint(int* a, int len) {
  for(int i = 0; i < len; i++)
    printf("%d, ", a[i]);
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
    //testExtract();
    int w = 16;

    int n = 7;
    int d = 6; // 2k-2 = 4
    int k = 3;
    int a = d - k + 1; // 4
    int B = k*a; // 12

    int* M = make_m_MSR(k, d);
    mprint(M, d, a, w);

    int* psi = talloc(int, n * d);
    int lambdas[n];
    make_psi_MSR_ext(psi, lambdas, n,k ,d, w); // n x d
    int* G = generator_for_PMC(psi, M, n, d, a, B);

    char* msg = "Regenerating codes are distributed storage codes\n\
that optimally trade bandwidth for storage...!";


    int msg_size = strlen(msg) + 1;
    printf("MESSAGE: \n%s\n", msg);
    printf("msg_size: %d\n", msg_size);

    char** all_dst_ptrs = talloc(char*, n*a);
    int chunk_size = msg_size / B; //bytes stored in all stripes for a particular coded-symbol.
    for (int i = 0; i < n*a; i++)
        all_dst_ptrs[i] = talloc(char, chunk_size);
    printf("chunk_size: %d\n", chunk_size);


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
    printf("Node 0 stores %d stripes of alpha=%d 2-byte-symbols each.\n", num_stripes, a);

    //cprint(node_ptrs, a, chunk_size);
    //printf("\n");


    /* Test DC. */
    int nodes[] = {0,1,3};
    char** chunk_ptrs_recv = talloc(char*, k*a);
    for (int i = 0; i < k; i++)
        for (int j = 0; j < a; j++)
            chunk_ptrs_recv[i*a + j] = all_dst_ptrs[nodes[i]*a + j];

    GeneratorDAG* dag = make_MSR_decoder_dag(psi, lambdas, M, nodes, k, d, w);
    print_dag(dag);

    char* dst = talloc(char, msg_size);
    apply_data_collect_dag_striped(dag, chunk_ptrs_recv, dst, B, msg_size, w);

    printf("decoded msg: \n");
    printf("%s\n", dst);
}
