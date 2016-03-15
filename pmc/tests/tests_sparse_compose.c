#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "jerasure.h"

#include "product_matrix.h"
#include "MBR_MSR.h"

#define talloc(type, num) (type *) malloc(sizeof(type)*(num))

#define mprint jerasure_print_matrix

void aprint(int* a, int len) {
  for(int i = 0; i < len; i++)
    printf("%d, ", a[i]);
  printf("\n");
}

int testLR_mult() {
  int w = 8;

  int A[12] =
   {1,2,3,
    4,5,6,
    7,8,9,
    10,11,12};

  // multiply A*X, where
  // A : 4 x 3
  // X : 3 x 2
  int n = 4;
  int k = 3;
  int m = 2;
  int* GL = left_mult_G(A, n, k, m);

  printf("G (AX):\n");
  mprint(GL, n*m, k*m, w);

  // mutiply X*A, where
  // X : 2 x 4
  // A : 4 x 3
  n = 2;
  k = 4;
  m = 3;
  int* GR = right_mult_G(A, n, k, m);
  printf("G (XA):\n");
  mprint(GR, n*m, n*k, w);

  free(GL);
  free(GR);
  return 0;
}
void test_single() {
    int w = 8;

    int n = 5;
    int d = 4;
    int k = 3;
    int a = d - k + 1; // 2
    int B = k*a;

    int* M = make_m_MSR(k, d);
    printf("M:\n");
    mprint(M, d, a, w);

    int* psi = talloc(int, n * d);
    int lambdas[n];
    make_psi_MSR_ext(psi, lambdas, n,k ,d, w); // n x d
    printf("psi:\n");
    mprint(psi, n, d, w);
    printf("lambdas: \n");
    aprint(lambdas, n);

    int* G = generator_for_PMC(psi, M, n, d, a, B);

    printf("G: \n");
    mprint(G, n*a, B, w);


    /* Encode msg to n nodes. */
    int msg[] = {0, 1, 2, 3, 4, 5};
    int* y = encode(msg, G, n, a, B, w); // (na)-vector
    printf("message: ");
    aprint(msg, B);
    printf("encoded: ");
    aprint(y, n*a);

    /* Test DC. */
    int nodes[] = {0,1,3};
    int* msg_recv = talloc(int, k*a);
    int* m_ptr = msg_recv;
    for(int i = 0; i < k; i++) {
        memcpy(m_ptr, y + a*nodes[i], sizeof(int) * a);
        m_ptr += a;
    }
    printf("recieved coded msg: ");
    aprint(msg_recv, B);

    GeneratorPlan* plan = make_MSR_2k2_decoder_plan(psi, lambdas, M, nodes, k, w);
    //printf("D0:\n");
    //mprint(plan->D[0], plan->dim[0], plan->dim[1], w);

    int* msg_dc = apply_plan_single(plan, msg_recv, w);
    printf("decoded msg: \n");
    aprint(msg_dc, B);
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
    int w = 16;

    int n = 5;
    int d = 4;
    int k = 3;
    int a = d - k + 1; // 2
    int B = k*a; // 6

    int* M = make_m_MSR(k, d);
    int* psi = talloc(int, n * d);
    int lambdas[n];
    make_psi_MSR_ext(psi, lambdas, n,k ,d, w); // n x d
    int* G = generator_for_PMC(psi, M, n, d, a, B);

    char* msg = "hello, world. striping test. last sentence now.";
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

    GeneratorPlan* plan = make_MSR_2k2_decoder_plan(psi, lambdas, M, nodes, k, w);

/*
 *    char** decoded_ptrs = talloc(char*, B);
 *    for (int i = 0; i < B; i++)
 *        decoded_ptrs[i] = talloc(char, chunk_size);
 *
 *    apply_plan_striped(plan, chunk_ptrs_recv, decoded_ptrs, chunk_size, w);
 *    printf("decoded chunks: \n");
 *    cprint(decoded_ptrs, B, chunk_size);
 */

    char* dst = talloc(char, msg_size);
    apply_data_collect_plan_striped(plan, chunk_ptrs_recv, dst, B, msg_size, w);
    printf("decoded msg: \n");
    printf("%s\n", dst);
}
