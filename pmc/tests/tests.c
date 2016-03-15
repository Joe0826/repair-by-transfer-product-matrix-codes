#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../jerasure-kmg/include/galois.h"
#include "reed_sol.h"
#include "jerasure.h"

#include "matrix.h"
#include "product_matrix.h"
#include "MBR_MSR.h"

#define talloc(type, num) (type *) malloc(sizeof(type)*(num))

#define mprint jerasure_print_matrix

void testRREF()
{
  int w = 16;

  int *a = reed_sol_vandermonde_coding_matrix(5,3,w);
  rref(a, 3, 5, w);
  printf("rref:\n");
  mprint(a, 3, 5, w);
  free(a);
}

void testRREF_aug()
{
  int w = 16;
  int A[16] =
   {1,0,1,0,
    0,1,1,1,
    1,1,1,1,
    1,0,1,1};
  rref_aug(A, 4, 4, 3, w);
  mprint(A, 4, 4, w);
}

void testCREF()
{
  int w = 16;

  int *a = reed_sol_vandermonde_coding_matrix(3,5,w);
  cref(a, 5, 3, w);
  printf("rref:\n");
  mprint(a, 5, 3, w);
  free(a);
}

void aprint(int* a, int len) {
  for(int i = 0; i < len; i++)
    printf("%d, ", a[i]);
  printf("\n");
}

int main(int argc, char **argv)
{
  int w = 16;

  const int n = 5;
  const int d = 4;
  const int a = 4;
  const int k = 3; //for M below:
  int M[16] =
    {1,2,3,7,
    2,4,5,8,
    3,5,6,9,
    7,8,9,0};
  int B = 9;

  int* psi = reed_sol_vandermonde_coding_matrix(d,n,w); // n x d
  mprint(psi, n, d, w);

  int* G = generator_for_PMC(psi, M, n, d, a, B);

  //printf("G:\n");
  //mprint(G, n*a, B, w);

  make_systematic(G, n, a, B, w);

  printf("G (systmatic):\n");
  mprint(G, n*a, B, w);


  /* Compute DC matrix. */
  int nodes[] = {0, 3, 4};
  int* D = make_data_collect_matrix(G, k, a, B, nodes, w);

  printf("D: \n");
  mprint(D, B, k*a, w);

  /* Encode msg to n nodes. */
  int msg[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  int* y = encode(msg, G, n, a, B, w); // (na)-vector
  printf("message: ");
  aprint(msg, B);

  /* Test DC. */
  int* yDC = talloc(int, k*a);
  for(int i = 0; i < k; i++) {
    memcpy(yDC + i*a, y + nodes[i]*a, sizeof(int)*a);
  }
  int* msg_dc = data_collect(D, yDC, k, a, B, w);
  printf("collected message: ");
  aprint(msg_dc, B);


  /* Test repair. */
  int f = 0;

  int* u_f = talloc(int, a);
  memcpy(u_f, psi + f*d, sizeof(int)*a); // only works for MBR

  int helper_nodes[d];
  int* yR = talloc(int, d);
  for(int i = 0, h = 0; h < d; i++) {
    if (i != f) {
      helper_nodes[h] = i;

      int* sym = helper_symbol(u_f, y + a*i, a, w);
      yR[h] = *sym;
      h++;
    }
  }


  int* A = make_repair_matrix(G, d, a, B, f, u_f, helper_nodes, w);
  printf("A: \n");
  mprint(A, a, d, w);
  printf("yR: ");
  aprint(yR, d);
  int* a_repaired = repair(A, yR, d, a, w);

  printf("stored data in node %d: \n", f);
  aprint(y + a*f, a);
  printf("repaired data: \n");
  aprint(a_repaired, a);


  free(psi);
  free(G);

  //testRREF();
  //testCREF();
  return 0;
}
