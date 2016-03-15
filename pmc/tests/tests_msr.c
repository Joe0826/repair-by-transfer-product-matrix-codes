#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../jerasure-kmg/include/galois.h"
#include "reed_sol.h"
#include "jerasure.h"

#include "product_matrix.h"
#include "MBR_MSR.h"

#define talloc(type, num) (type *) malloc(sizeof(type)*(num))

#define mprint jerasure_print_matrix

void testDC()
{
  //TODO
}

void aprint(int* a, int len) {
  for(int i = 0; i < len; i++)
    printf("%d, ", a[i]);
  printf("\n");
}

void test0() {
  int w = 16;
  int k = 3;
  int d = 6;
  int n = 7;
  int* message_matrix = make_m_MSR(k, d);
  printf("matrix M:\n");
  mprint(message_matrix, d, d - k + 1, w);
  free(message_matrix);


  int* psi_matrix = make_psi_MSR(n, k, d, w);
  printf("matrix psi:\n");
  mprint(psi_matrix, n, d, w);
  free(psi_matrix);

  int* mu;
  for (int i = 0; i < n;i++) {
    printf("row %d:\n", i);
    mu = make_mu_MSR(i, k, d, w);
	mprint(mu, 1, d, w);
	free(mu);
  }
}

void test1() {
  int w = 16;

  int n = 5;
  int d = 4;
  int k = 2;
  int a = 3;

  int* M = make_m_MSR(k, d);
  printf("M:\n");
  mprint(M, d, d - k + 1, w);
  int B = 6;

  int* psi = make_psi_MSR(n,k ,d, w); // n x d
  printf("psi:\n");
  mprint(psi, n, d, w);

  int* G = generator_for_PMC(psi, M, n, d, a, B);

  make_systematic(G, n, a, B, w);

  printf("G (systmatic):\n");
  mprint(G, n*a, B, w);


  /* Compute DC matrix. */
  int nodes[] = {0, 3};
  int* D = make_data_collect_matrix(G, k, a, B, nodes, w);

  printf("D: \n");
  mprint(D, B, k*a, w);

  /* Encode msg to n nodes. */
  int msg[] = {0, 1, 2, 3, 4, 5};
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

  //Repair
  int helper_nodes[d];
  int* yR = talloc(int, d);
  for (int f = 0; f < n; f++) {
    printf("Node %d failed!\n", f);
    int* u_f = make_mu_MSR(f, k, d, w);
	for (int i = 0, h = 0; h < d; i++) {
      if (i != f) {
	    helper_nodes[h] = i;
		int* sym = helper_symbol(u_f, y + a * i, a, w);
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
	free(u_f);
  }
  free(yR);
  free(M);
  free(psi);
  free(G);

}

void test2() {
  int w = 16;

  int n = 10;
  int d = 8;
  int a = 6;
  int k = 3; //for M below:
  int* M = make_m_MSR(k, d);
  int B = 18;


  printf("M:\n");
  mprint(M, d, d - k + 1, w);
  int* psi = make_psi_MSR(n, k, d, w); // n x d
  printf("psi:\n");
  mprint(psi, n, d, w);

  int* G = generator_for_PMC(psi, M, n, d, a, B);

  //printf("G:\n");
  //mprint(G, n*a, B, w);

  make_systematic(G, n, a, B, w);

  printf("G (systmatic):\n");
  mprint(G, n*a, B, w);


  /* Compute DC matrix. */
  int nodes[] = {0, 4, 6};
  int* D = make_data_collect_matrix(G, k, a, B, nodes, w);

  printf("D: \n");
  mprint(D, B, k*a, w);

  /* Encode msg to n nodes. */
  int msg[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17};
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


  free(M);
  free(psi);
  free(G);


}


int main(int argc, char **argv) {
  //test0();
  test1();
  test2();
  return 0;
}

