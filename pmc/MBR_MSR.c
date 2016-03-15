#include <string.h>
#include <stdio.h>
#include "../jerasure-kmg/include/galois.h"
#include "../jerasure-kmg/include/jerasure.h"
#include "matrix.h"
#include "product_matrix.h"
#include "dag_helper.h"
#include "MBR_MSR.h"

#define talloc(type, num) (type *) malloc(sizeof(type)*(num))

GeneratorDAG* make_MBR_decoder_dag(int* psi, int* M,
        int* nodes, int k, int d, int w) {
    int a = d;
    int B = k*d - k*(k-1)/2;

    int* phi = talloc(int, k*k);
    int* delta = talloc(int, k*(d-k));
    int *phi_r = phi, *delta_r = delta;
    for(int i = 0; i < k; i++) {
        memcpy(phi_r, psi + d*nodes[i], sizeof(int)*k);
        memcpy(delta_r, psi + k + d*nodes[i], sizeof(int)*(d-k));
        phi_r += k;
        delta_r += (d-k);
    }

    // extract CR (right half of C)
    GeneratorDAG* CR = extract_from_dag(NULL, k, a, 0, k, k, d-k, "C_{DC}^R");

    // multiply phi^-1 * CR to get T
    int* phi_inv = talloc(int, k*k);
    invert_matrix(phi, phi_inv, k, w);
    GeneratorDAG* T = left_mult_dag(phi_inv, CR, k, k, d-k, "T");
    GeneratorDAG* Tt = transpose_dag(T, k, d-k);

    GeneratorDAG* CL = extract_from_dag(NULL, k, a, 0, 0, k, k, "C_{DC}^L");

    GeneratorDAG* deltaTt = left_mult_dag(delta, Tt, k, d-k, k, "deltaTt");
    GeneratorDAG* phiS = sub_dag(CL, deltaTt, k, k, "phiS");
    GeneratorDAG* S = left_mult_dag(phi_inv, phiS, k, k, k, "S");



    // combine decoded submatrices, unwrap final message
    GeneratorDAG* L = vstack_dag(S, Tt, k, d-k, k, "M_L"); // left block of message mtx
    GeneratorDAG* R = vstack_dag(
            T, noop_dag((d-k)*(d-k)),
            k, d-k, d-k, "M_R");

    GeneratorDAG* M_dag = hstack_dag(L, R, k, d-k, d, "M");

    // finally, unwrap M into the message-vector
    return msg_unwrap_dag(M_dag, M, d, a, B);
}

GeneratorDAG* make_MSR_decoder_dag(int* psi, int* lambdas, int* M,
        int* nodes, int k, int d, int w) {
    if (d == 2*k - 2) {
        return make_MSR_2k2_decoder_dag(NULL, psi, lambdas, M, nodes, k, w);
    }

    int a = d - k + 1;
    int B = k*a;
    int u = d - (2*k-2);

    // NOTE: using unified-MSR, where psi = [lambda phi | phi | delta]
    int* lambdas_full = lambdas;
    lambdas = talloc(int, k);
    int* phi = talloc(int, k*(k-1));
    int* phi_hat = talloc(int, k*k);
    int* delta = talloc(int, k*u);
    int* delta_hat = talloc(int, k*(u-1));
    int *phi_r = phi, *phi_hat_r = phi_hat, *delta_hat_r = delta_hat, *delta_r = delta;
    for(int i = 0; i < k; i++) {
        lambdas[i] = lambdas_full[nodes[i]];
        memcpy(phi_r, psi + d*nodes[i] + (k-1), sizeof(int)*(k-1));
        memcpy(phi_hat_r, psi + d*nodes[i] + (k-1), sizeof(int)*k);
        memcpy(delta_r, psi + d*nodes[i] + (2*k - 2), sizeof(int)*u);
        memcpy(delta_hat_r, psi + d*nodes[i] + (2*k - 1), sizeof(int)*(u-1));
        phi_r += (k-1); phi_hat_r += k; delta_hat_r += u-1; delta_r += u;
    }

    // extract CR
    GeneratorDAG* CR = extract_from_dag(NULL, k, a, 0, k-1, k, u, "C_{DC}^R");

    // multiply phi^-1 * CR
    int* phi_hat_inv = talloc(int, k*k);
    jerasure_invert_matrix(phi_hat, phi_hat_inv, k, w);

    GeneratorDAG* pInvCR = left_mult_dag(phi_hat_inv, CR, k, k, u, "phi_hat^{-1} * CR");

    // extract z1, compute phi^{-1}*delta*z1
    GeneratorDAG* Z1 = extract_from_dag(pInvCR, k, u, k-1, 1, 1, u-1, "Z1");

    int* pInvDelta = talloc(int, k*(u-1));
    matrix_multiply(pInvDelta, phi_hat_inv, delta_hat, k, k, k, u-1, w);
    GeneratorDAG* pdZ1 = left_mult_dag(pInvDelta, Z1, k, u-1, 1, "phi_hat^{-1}*delta_hat* z1");


    // recover T_hat

    // first column of phi^{-1}*CR
    GeneratorDAG* pInvCR_1 = extract_from_dag(pInvCR, k, u, 0, 0, k, 1, "pInvCR_1");

    // first column of T_hat
    GeneratorDAG* T_hat_1 = add_dag(pInvCR_1, pdZ1, k, 1, "T_hat_1");

    GeneratorDAG* T_hat_tail = extract_from_dag(pInvCR, k, u, 0, 1, k, u-1, "T_hat_tail");

    GeneratorDAG* T_hat = hstack_dag(T_hat_1, T_hat_tail, 1, u-1, k, "T_hat");


    // extract CL, T, find C_msr = lambda*Phi*S1 + Phi*S2
    GeneratorDAG* CL = extract_from_dag(NULL, k, a, 0, 0, k, k-1, "C_{DC}^L");

    GeneratorDAG* T = extract_from_dag(T_hat, k, u, 0, 0, k-1, u, "T");
    GeneratorDAG* Tt = transpose_dag(T, k-1, u);
    GeneratorDAG* dTt = left_mult_dag(delta, Tt, k, u, k-1, "delta*T^T");

    GeneratorDAG* C_msr = add_dag(CL, dTt, k, k-1, "C_msr");


    // legacy MSR plan (hack) TODO: unhack
    int* M_fake = calloc(2*(k-1)*(k-1), sizeof(int));

    int* psi_msr = talloc(int, (2*k-2)*k);
    int* nodes_msr = talloc(int, k);
    for(int i = 0; i < k; i++) {
        nodes_msr[i] = i;
        memcpy(psi_msr + i*(2*k-2), psi + nodes[i]*d, sizeof(int)*(2*k-2));
    }

    GeneratorDAG* msr_dag = make_MSR_2k2_decoder_dag(C_msr, psi_msr, lambdas, M_fake, nodes_msr, k, w);
    msr_dag = msr_dag->deps[0]; // remove the last stage, ie don't unwrap [S1 \\ S2] (hack)

    // combine decoded submatrices, unwrap final message

    // construct L, R, --> M matrix (vhstacks)
    // then write the inverse-extract fund

    GeneratorDAG* L = vstack_dag(msr_dag, Tt, 2*k-2, u, k-1, "M_L"); // left block of message mtx

    GeneratorDAG* Z_hat = hstack_dag(
            Z1, noop_dag((u-1)*(u-1)),
            1, u-1,
            u-1, "Z_hat");
    GeneratorDAG* TZ = vstack_dag(T_hat, Z_hat, k, u-1, u, "TZ");
    GeneratorDAG* R = vstack_dag(
            noop_dag((k-1)*u), TZ,
            k-1, (k-1)+u,
            u, "M_R");

    GeneratorDAG* M_dag = hstack_dag(L, R, k-1, u, d, "M");

    // finally, unwrap M into the message-vector
    return msg_unwrap_dag(M_dag, M, d, a, B);
}

GeneratorDAG* make_MSR_2k2_decoder_dag(GeneratorDAG* src, int* psi, int* lambdas, int* M,
        int* nodes, int k, int w) {

    int a = k-1;
    int d = 2*a;
    int B = k*a;

    // NOTE: using unified-MSR, where psi = [lambda phi | phi]
    int* lambs = talloc(int, k);
    int* phi = talloc(int, k*a);
    int* phi_r = phi;
    for(int i = 0; i < k; i++) {
        lambs[i] = lambdas[nodes[i]];
        memcpy(phi_r, psi + d*nodes[i] + a, sizeof(int)*a);
        phi_r += a;
    }
    lambdas = lambs;
    int* phiT = talloc(int, a*k);
    transpose(phi, phiT, k, a);

    // takes the encoded message (psi M : k x a)
    // to X = Lambda P + Q : (k x k)
    int* D1 = right_mult_G(phiT, k, a, k);
    GeneratorDAG* stage1 = make_transform_dag(src, D1, k*k, "2k2_stage1");

    //define P' = P, discarding diagonal entries. k x (k-1)

    // we construct D2 to take us from vec(X) --> vec( vstack(P', Q') )
    // conveniently, this is the same as vstack( vec(P'), vec(Q') )

    // to construct P' and Q' : 2k x (k-1)   from X : k x k
    int* D2 = talloc(int, (2*k)*(k-1) * (k*k));
    memset(D2, 0, sizeof(int) * (2*k)*(k-1) * (k*k));
    // (i, j): indices into X
    for (int i = 0; i < k; i++) {
        int ip = i; // indices into P'
        int jp = 0;
        for (int j = 0; j < k; j++) {
            if (i==j) continue;


            // P_{ij} = (1 / (l_i - l_j)) x_{ij} - (1 / (l_i - l_j)) x_{ji}
            int dLambda = lambdas[i] ^ lambdas[j]; //l_i - l_j
            int eP = ((k-1)*ip + jp) * (k*k); // row corrosponding to P_{ij} : (k-1)*ip + jp'th row
            // so eP'th entry in vec( vstack(P', Q') ) = P_{ij}
            D2[eP + (k*i + j)] = D2[eP + (k*j + i)] = galois_single_divide(1, dLambda, w);


            // Q_{ij} = (-l_j / (l_i - l_j)) x_{ij} + (l_i / (l_i - l_j)) x_{ji}
            int eQ = eP + k*(k-1) * (k*k); // we offset by len(vec(P')) rows.
            // so eQ'th entry in vec( vstack(P', Q') ) = Q_{ij}
            D2[eQ + (k*i + j)] = galois_single_divide(lambdas[j], dLambda, w); // x_ij factor
            D2[eQ + (k*j + i)] = galois_single_divide(lambdas[i], dLambda, w); // x_ji factor

            jp++;
        }
    }
    GeneratorDAG* stage2 = make_transform_dag(stage1, D2, 2*k*(k-1), "2k2_stage2_PQ");

    // total D for taking [P' \\ Q'] --> [psi{a} * S1 \\ psi{a} * S2]
    int D3_cols = 2*k*(k-1);
    int* D3 = talloc(int, (d*a) * D3_cols);
    memset(D3, 0, sizeof(int)*(d*a)*D3_cols);

    int* phi_a = talloc(int, a*a); // phi, with the i'th row removed
    int* phi_aT = talloc(int, a*a);
    int* phi_aT_Inv = talloc(int, a*a);
    int* D_i = talloc(int, a*a); // D for computing row i of [psi * S_1] (1 x a) from row i of P' (1 x k-1=a)
    int* D3_r = D3; // beginning of block for vec( row i of psi * S_1 )

    for (int i = 0; i < a; i++) {
        // setup phi_a
        int* phi_r = phi;
        int* phi_ar = phi_a;
        for (int r = 0; r < a+1; r++) {
            if (r == i) {
                phi_r += a;
                continue;
            }
            memcpy(phi_ar, phi_r, sizeof(int)*a);
            phi_r += a;
            phi_ar += a;
        }

        transpose(phi_a, phi_aT, a, a);

        jerasure_invert_matrix(phi_aT, phi_aT_Inv, a, w);

        // compute D_i
        right_mult_G_into(D_i, phi_aT_Inv, 1, a, a);

        // place D_i into D3
        // D3 = [  ...  ]
        //      [D_i   0]
        //      [  ...  ]
        //      [0   D_i]
        //printf("D_i:\n");
        //jerasure_print_matrix(D_i, a, a, w);

        submat_cpy(D_i, a, a, a, D3_r + i*a, D3_cols); // for P'/S1
        submat_cpy(D_i, a, a, a, D3_r + i*a + a*a*D3_cols + k*(k-1), D3_cols); // for Q'/S2

        D3_r += a*D3_cols;
    }
    GeneratorDAG* stage3 = make_transform_dag(stage2, D3, d*a, "2k2_stage3_PsiS");
    stage3 = merge_top2_transforms(stage3, w);

    // to take [psi{a} * S1 \\ psi{a} * S2] --> [S1 \\ S2]
    int* D4 = talloc(int, (d*a) * (d*a));
    memset(D4, 0, sizeof(int)*(d*a) * (d*a));

    // setup phi_a = first alpha rows of phi
    phi_r = phi;
    int* phi_ar = phi_a;
    for (int r = 0; r < a; r++) {
        memcpy(phi_ar, phi_r, sizeof(int)*a);
        phi_r += a;
        phi_ar += a;
    }
    int* phi_a_Inv = phi_aT_Inv; //re-use memory //talloc(int, a*a);
    jerasure_invert_matrix(phi_a, phi_a_Inv, a, w);

    //printf("phi_a_Inv:\n");
    //jerasure_print_matrix(phi_a_Inv, a, a, w);

    // compute inverse, copy into D4 for S1 and S2
    // D4 = [D4_1   0]
    //      [0   D4_1]
    int* D4_1 = left_mult_G(phi_a_Inv, a, a, a);
    submat_cpy(D4_1, a*a, a*a, a*a, D4, d*a);
    submat_cpy(D4_1, a*a, a*a, a*a, D4 + (a*a)*(d*a) + a*a, d*a);
    GeneratorDAG* stage4 = make_transform_dag(stage3, D4, d*a, "2k2_stage4_S");

    // finally, unpack [S1 \\ S2] into the original message vector
    GeneratorDAG* final = msg_unwrap_dag(stage4, M, d, a, B);
    return final;
}

GeneratorPlan* make_MSR_2k2_decoder_plan(int* psi, int* lambdas, int* M,
        int* nodes, int k, int w) {

    int a = k-1;
    int d = 2*a;
    int B = k*a;

    GeneratorPlan* plan = (GeneratorPlan*) malloc(sizeof(GeneratorPlan));
    int numD = 0;


    // NOTE: using unified-MSR, where psi = [lambda phi | phi]
    int* lambs = talloc(int, k);
    int* phi = talloc(int, k*a);
    int* phi_r = phi;
    for(int i = 0; i < k; i++) {
        lambs[i] = lambdas[nodes[i]];
        memcpy(phi_r, psi + d*nodes[i] + a, sizeof(int)*a);
        phi_r += a;
    }
    lambdas = lambs;
    int* phiT = talloc(int, a*k);
    transpose(phi, phiT, k, a);

    // takes the encoded message (psi M : k x a)
    // to X = Lambda P + Q : (k x k)
    int* D1 = right_mult_G(phiT, k, a, k);
    plan->dim[numD] = k*a;    // input size
    plan->dim[numD+1] = k*k;  // output size
    plan->D[numD++] = D1;

    //define P' = P, discarding diagonal entries. k x (k-1)

    // we construct D2 to take us from vec(X) --> vec( vstack(P', Q') )
    // conveniently, this is the same as vstack( vec(P'), vec(Q') )

    // to construct P' and Q' : 2k x (k-1)   from X : k x k
    int* D2 = talloc(int, (2*k)*(k-1) * (k*k));
    memset(D2, 0, sizeof(int) * (2*k)*(k-1) * (k*k));
    // (i, j): indices into X
    for (int i = 0; i < k; i++) {
        int ip = i; // indices into P'
        int jp = 0;
        for (int j = 0; j < k; j++) {
            if (i==j) continue;


            // P_{ij} = (1 / (l_i - l_j)) x_{ij} - (1 / (l_i - l_j)) x_{ji}
            int dLambda = lambdas[i] ^ lambdas[j]; //l_i - l_j
            int eP = ((k-1)*ip + jp) * (k*k); // row corrosponding to P_{ij} : (k-1)*ip + jp'th row
            // so eP'th entry in vec( vstack(P', Q') ) = P_{ij}
            D2[eP + (k*i + j)] = D2[eP + (k*j + i)] = galois_single_divide(1, dLambda, w);


            // Q_{ij} = (-l_j / (l_i - l_j)) x_{ij} + (l_i / (l_i - l_j)) x_{ji}
            int eQ = eP + k*(k-1) * (k*k); // we offset by len(vec(P')) rows.
            // so eQ'th entry in vec( vstack(P', Q') ) = Q_{ij}
            D2[eQ + (k*i + j)] = galois_single_divide(lambdas[j], dLambda, w); // x_ij factor
            D2[eQ + (k*j + i)] = galois_single_divide(lambdas[i], dLambda, w); // x_ji factor

            jp++;
        }
    }
    plan->D[numD++] = D2;
    plan->dim[numD] = 2*k*(k-1);

    // total D for taking [P' \\ Q'] --> [psi{a} * S1 \\ psi{a} * S2]
    int D3_cols = 2*k*(k-1);
    int* D3 = talloc(int, (d*a) * D3_cols);
    memset(D3, 0, sizeof(int)*(d*a)*D3_cols);

    int* phi_a = talloc(int, a*a); // phi, with the i'th row removed
    int* phi_aT = talloc(int, a*a);
    int* phi_aT_Inv = talloc(int, a*a);
    int* D_i = talloc(int, a*a); // D for computing row i of [psi * S_1] (1 x a) from row i of P' (1 x k-1=a)
    int* D3_r = D3; // beginning of block for vec( row i of psi * S_1 )

    for (int i = 0; i < a; i++) {
        // setup phi_a
        int* phi_r = phi;
        int* phi_ar = phi_a;
        for (int r = 0; r < a+1; r++) {
            if (r == i) {
                phi_r += a;
                continue;
            }
            memcpy(phi_ar, phi_r, sizeof(int)*a);
            phi_r += a;
            phi_ar += a;
        }

        transpose(phi_a, phi_aT, a, a);

        jerasure_invert_matrix(phi_aT, phi_aT_Inv, a, w);

        // compute D_i
        right_mult_G_into(D_i, phi_aT_Inv, 1, a, a);

        // place D_i into D3
        // D3 = [  ...  ]
        //      [D_i   0]
        //      [  ...  ]
        //      [0   D_i]
        //printf("D_i:\n");
        //jerasure_print_matrix(D_i, a, a, w);

        submat_cpy(D_i, a, a, a, D3_r + i*a, D3_cols); // for P'/S1
        submat_cpy(D_i, a, a, a, D3_r + i*a + a*a*D3_cols + k*(k-1), D3_cols); // for Q'/S2

        D3_r += a*D3_cols;
    }
    plan->D[numD++] = D3;
    plan->dim[numD] = d*a;


    // to take [psi{a} * S1 \\ psi{a} * S2] --> [S1 \\ S2]
    int* D4 = talloc(int, (d*a) * (d*a));
    memset(D4, 0, sizeof(int)*(d*a) * (d*a));

    // setup phi_a = first alpha rows of phi
    phi_r = phi;
    int* phi_ar = phi_a;
    for (int r = 0; r < a; r++) {
        memcpy(phi_ar, phi_r, sizeof(int)*a);
        phi_r += a;
        phi_ar += a;
    }
    int* phi_a_Inv = phi_aT_Inv; //re-use memory //talloc(int, a*a);
    jerasure_invert_matrix(phi_a, phi_a_Inv, a, w);

    //printf("phi_a_Inv:\n");
    //jerasure_print_matrix(phi_a_Inv, a, a, w);

    // compute inverse, copy into D4 for S1 and S2
    // D4 = [D4_1   0]
    //      [0   D4_1]
    int* D4_1 = left_mult_G(phi_a_Inv, a, a, a);
    submat_cpy(D4_1, a*a, a*a, a*a, D4, d*a);
    submat_cpy(D4_1, a*a, a*a, a*a, D4 + (a*a)*(d*a) + a*a, d*a);

    plan->D[numD++] = D4;
    plan->dim[numD] = d*a;

    // finally, unpack [S1 \\ S2] into the original message vector
    int* D5 = talloc(int, B * (d*a));
    memset(D5, 0, sizeof(int)*B*(d*a));
    char* hit = (char*) malloc(sizeof(char) * B);
    memset(hit, 0, sizeof(char)*B);
    for (int i = 0; i < d*a; i++) {
        int e = M[i] - 1;
        if (e >= 0 && !hit[e]) {
            D5[e * (d*a) + i] = 1;
            hit[e] = 1;
        }
    }
    plan->D[numD++] = D5;
    plan->dim[numD] = B;

    // TODO:
    // - collapse unpacking and last decoder step
    // - free

    plan->D[numD] = NULL;
    plan->dim[numD+1] = 0;
    return plan;
}

/*
  Generate matrix M
*/
int* make_m_MBR(int k, int d) {
  int elem = 1;
  int* M = talloc(int, d * d);
  memset(M, 0, sizeof(int)*(d*d));
  // Fill in the upper triangular region of S
  for (int i = 0; i < k; i++) {
    for (int j = i; j < k; j++) {
      M[i * d + j] = elem;
      M[j * d + i] = elem;
      elem++;
    }
  }
  // Fill in T
  for (int i = 0; i < k; i++) {
    for (int j = k; j < d; j++) {
      M[i * d + j] = elem;
      elem++;
    }
  }
  // Fill in T_t
  for (int i = k; i < d; i++) {
    for (int j = 0; j < k; j++) {
      M[i * d + j] = M[j * d + i];
    }
  }
  // Fill in 0
  for (int i = k; i < d; i++) {
    for (int j = k; j < d; j++) {
      M[i * d + j] = 0;
    }
  }
  return M;
}

int* make_psi_MBR(int n, int k, int d, int w) {
  int* psi = talloc(int, n * d);

  for (int i = 0; i < k; i++) {
    for (int j = 0; j < d; j++) {
      psi[i * d + j] = 0;
    }
    psi[i * d + i] = 1;
  }
  for (int i = k; i < n; i++) {
    int x_i = i + d - k;
    for (int j = 0; j < d; j++) {
      int y_j = j;
      psi[i * d + j] = galois_single_divide(1, x_i ^ y_j, w); // x_i - y_j
    }
  }

  return psi;
}

void make_mu_MBR_inplace(int* mu, const int f, const int k, const int d, const int w) {
  if (f < k) {
    for (int j = 0; j < d; j++) {
      mu[j] = 0;
    }
    mu[f] = 1;
  } else {
    int x_i = f + d - k;
    for (int j = 0; j < d; j++) {
      int y_j = j;
      mu[j] = galois_single_divide(1, x_i ^ y_j, w); // x_i - y_j
    }
  }
}

/*
  Generate mu_f (the helper-vector for node f failing).
*/
int* make_mu_MBR(int f, int k, int d, int w) {
  int* mu = talloc(int, d);
  make_mu_MBR_inplace(mu, f, k, d, w);
  return mu;
}
/*
Generate d-by-(d - k + 1) matrix
*/
int* make_m_MSR(int k, int d) {
  int elem = 1;
  int num_col = d - k + 1;
  int* M = talloc(int, d * (d - k + 1));
  memset(M, 0, sizeof(int) * d * (d - k + 1));
  // Fill in the symmetric matrix S1
  for (int i = 0; i < k - 1; i++) {
    for (int j = i; j < k - 1; j++) {
      M[i * num_col + j] = elem;
      M[j * num_col + i] = elem;
      elem++;
    }
  }
  // Fill in the symmetric matrix S2
  for (int i = k - 1; i < 2 * (k - 1); i++) {
    for (int j = i - (k - 1); j < k - 1; j++) {
      M[i * num_col + j] = elem;
      M[(k - 1 + j) * num_col + i - (k - 1)] = elem;
      elem++;
    }
  }
  // Fill in T
  for (int i = k - 1; i < 2 * (k - 1); i++) {
    for (int j = k - 1; j < d - k + 1; j++) {
      M[i * num_col + j] = elem;
      elem++;
    }
  }
  // Fill in T_t
  for (int i = 2 * (k - 1); i < d; i++) {
    for (int j = 0; j < k - 1; j++) {
      M[i * num_col + j] = M[(k - 1 + j) * num_col + (i - (k - 1))];
    }
  }
  // Fill in Z
  // Fill in z_0 and z_1^t
  for (int j = k - 1; j < d - k + 1; j++) {
    M[2 * (k - 1) * num_col + j] = elem;
    elem++;
  }
  // Fill in z_1
  for (int i = 2 * (k - 1); i < d; i++) {
      M[i * num_col + k - 1] = M[2 * (k - 1) * num_col + i - ( k - 1)];
  }
  return M;
}

int* make_psi_MSR(int n, int k, int d, int w) {
    int* psi = talloc(int, n*d);
    int* lambdas = talloc(int, n);
    make_psi_MSR_ext(psi, lambdas, n, k, d, w);
    return psi;
}

void make_psi_MSR_ext(int* psi, int* lambdas, int n, int k, int d, int w) {
  // psi should be n x d

  for (int i = 0; i < n; i++) {
    // Fill in the (k-1)-th column 1
    psi[i * d + k - 1] = 1;
    for (int j = 1; j < 2 * k - 2; j++) {
      if (j % 2 == 0) {
        psi[i * d + k - 1 + j / 2] = galois_single_multiply(psi[i * d + j / 2 - 1], i + 1, w);
      } else {
        psi[i * d + (j - 1) / 2] = galois_single_multiply(psi[i * d + (j - 1) / 2 + k - 1], i + 1, w);
      }
    }
    lambdas[i] = psi[i * d];
  }
  if (d > 2*k - 2) {
      // Fill in Delta.
      for (int i = 0; i < n; i++) {
        int square = galois_single_multiply(i+1, i+1, w);
        psi[i * d + 2*k - 2] = galois_single_multiply(psi[i * d + 2*k - 3], square, w);
        for (int j = 2 * k - 1; j < d; j++) {
          psi[i * d + j] = galois_single_multiply(psi[i * d + j - 1], i + 1, w);
        }
      }
  }
}

void make_mu_MSR_inplace(int* mu, const int f, const int k, const int d, const int w) {
  mu[0] = 1;
  int square = galois_single_multiply(f+1, f+1, w);
  for (int i = 1; i < k - 1; i++) {
    mu[i] = galois_single_multiply(mu[i - 1], square, w);
  }
  if (d > 2*k - 2) {
      mu[k-1] = galois_single_multiply(mu[k-2], square, w);
      for (int i = k; i < d - k + 1; i++) {
        mu[i] = galois_single_multiply(mu[i - 1], f + 1, w);
      }
  }
}

/*
  Generate mu_f (the helper-vector for node f failing).
*/
int* make_mu_MSR(int f, int k, int d, int w) {
  int* mu = talloc(int, d - k + 1);
  make_mu_MSR_inplace(mu, f, k, d, w);
  return mu;
}
