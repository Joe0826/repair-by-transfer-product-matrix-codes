#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../jerasure-kmg/include/jerasure.h"

#include "test_common.h"
#include "matrix.h"
#include "product_matrix.h"
#include "dag_helper.h"
#include "MBR_MSR.h"
#include "rbt.h"

int main(int argc, char **argv)
{
    const int w = 16;

    const int n = 7;
    const int d = 6; // 2k-2 = 4
    const int k = 3;
    const int a = d - k + 1; // 4
    const int B = k*a; // 12

    int* M = make_m_MSR(k, d);
    //mprint(M, d, a, w);

    int* psi = talloc(int, n * d);
    int lambdas[n];
    make_psi_MSR_ext(psi, lambdas, n, k, d, w); // n x d
    int* G = generator_for_PMC(psi, M, n, d, a, B);

    char* msg = "Regenerating codes are distributed storage codes\n\
that optimally trade bandwidth for storage...!";


    int msg_size = strlen(msg) + 1;
    printf("MESSAGE: \n%s\n", msg);
    printf("msg_size: %d\n", msg_size);

    char** all_dst_ptrs = talloc(char*, n*a);
    int chunk_size = msg_size / B;  // bytes stored in all stripes for a particular coded-symbol.
    for (int i = 0; i < n*a; i++)
        all_dst_ptrs[i] = talloc(char, chunk_size);
    printf("chunk_size: %d\n", chunk_size);


    /* Prepare RBT (repair-by-transfer) data. */
    int* all_mus = talloc(int, n*a);
    for (int i = 0; i < n; i++) {
        make_mu_MSR_inplace(all_mus + i*a, i, k, d, w);
    }
    int* all_P;
    int* all_P_inv;
    int* nodes_P_idx;
    gen_rbt_data(&all_P, &all_P_inv, &nodes_P_idx, all_mus, n, a, w, RBT_CYC);

    /* Encode all nodes. */
    char** src_ptrs = init_data_ptrs(msg, B, msg_size);
    char** dst_i_ptr = all_dst_ptrs;
    int* G_i = G;
    for (int i = 0; i < n; i++) {
        // first normal encoding (regular generator mtx multiply)
        // the DAG equivalent of:
        // encode_striped(msg, G_i, dst_i_ptr, a, B, msg_size, w);
        GeneratorDAG* encode_dag = make_transform_dag(NULL, G_i, a, "pmc_encode");

        // then RBT-postcoding
        int* P = all_P + (a*a)*nodes_P_idx[i];
        GeneratorDAG* rbt_post = make_transform_dag(encode_dag, P, a, "rbt_postcode");

        // perform coding.
        apply_dag_striped(rbt_post, src_ptrs, dst_i_ptr, B, chunk_size, w);

        // move to next node
        dst_i_ptr += a;
        G_i += a*B;
    }


    int num_stripes = msg_size / (B * w/8);
    char** all_data_ptrs = all_dst_ptrs;
    printf("Node 0 stores %d stripes of alpha=%d 2-byte-symbols each.\n",
            num_stripes, a);

    /* Test Repair. */

    // Repair-By-Transfer
    // for fixing node 0 from nodes 1-6.
    // Nodes 3, 4, 5, 6 will simply transfer data (RBT),
    // Nodes 1, 2 will compute helper-chunks.
    int f_id = 0;

    char* helper_chunk_ptrs[d];
    int helper_nodes[] = {1, 2, 3, 4, 5, 6};
    int* orig_mu_f = all_mus + f_id * a;  // make_mu_MSR(f_id, k, d, w);

    // Generate (or transfer) helper-chunks.
    for (int i = 0; i < d; i++) {
        int h_id = helper_nodes[i];
        int idx = get_repair_chunk_idx(h_id, f_id, n, a, RBT_CYC);
        if (idx >= 0) {
            // We have the exact repair chunk, can just copy.
            helper_chunk_ptrs[i] = all_data_ptrs[h_id*a + idx];
        } else {
            // Need to compute the repair chunk.
            int* P_inv = all_P_inv + (a*a)*nodes_P_idx[h_id];
            int* mu_f = rbt_mu_f(orig_mu_f, P_inv, a, w);
            char** chunk_ptrs = all_data_ptrs + h_id*a; // stored chunks.
            char* helper_chunk = talloc(char, chunk_size);
            helper_symbol_striped(chunk_ptrs, mu_f, helper_chunk, a, B, msg_size, w);
            helper_chunk_ptrs[i] = helper_chunk;
        }
    }

    // Repair.
    int* Pf = all_P + (a*a)*f_id;
    int* A_orig = make_repair_matrix(G, d, a, B, f_id, orig_mu_f, helper_nodes, w);

    // First repair as usual, then rbt-postcode: A = Pf * A_orig
    int* A = talloc(int, a * d);
    matrix_multiply(A, Pf, A_orig, a, a, a, d, w);

    // Alloc repaired chunk
    char* repaired_data_ptrs[a];
    for (int i = 0; i < a; i++) {
        repaired_data_ptrs[i] = (char*) malloc(chunk_size);
    }
    // Compute repaired chunk
    repair_striped(helper_chunk_ptrs, A, repaired_data_ptrs,
            d, a, B, msg_size, w);

    // Check validity.
    for (int i = 0; i < a; i++)
        for (int b = 0; b < chunk_size; b++)
            if (repaired_data_ptrs[i][b] != all_data_ptrs[f_id*a + i][b]) {
                fprintf(stderr, "ERROR: repair invalid.\n");
                return 1;
            }
    printf("\nSUCCESS repairing failed node!\n\n");

    // Mess with memory, to try to break
    // any assumptions about uninitialized memory.
    memory_monkey();

    /* Test Data-Collect. */

    int dc_nodes[] = {0,1,3};
    char** chunk_ptrs_recv = talloc(char*, k*a);
    for (int i = 0; i < k; i++)
        for (int j = 0; j < a; j++)
            chunk_ptrs_recv[i*a + j] = all_data_ptrs[dc_nodes[i]*a + j];


    // Unshuffle RBT.
    GeneratorDAG* rbt_decode = rbt_decode_dc(all_P_inv, dc_nodes, k, a, w);
    print_dag(rbt_decode);

    // Normal dc-decoder.
    GeneratorDAG* dc_dag = make_MSR_decoder_dag(psi, lambdas, M, dc_nodes, k, d, w);

    // First rbt_decode, then normal data-collect.
    GeneratorDAG* full_dc_dag = concat_dag(rbt_decode, dc_dag);

    char* dst = talloc(char, msg_size);
    apply_data_collect_dag_striped(full_dc_dag, chunk_ptrs_recv, dst, k, a, B, msg_size, w);

    printf("\ndecoded msg: \n");
    printf("%s\n", dst);
}
