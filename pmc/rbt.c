#include "rbt.h"

#include <string.h>
#include "../jerasure-kmg/include/jerasure.h"
#include "matrix.h"
#include "dag_helper.h"

/*
 * Returns the id of the i'th child of given node_id,
 * accodring to the RbtType.
 */
int get_child_id(const int node_id, const int i, const int n,
        const RbtType rbt) {
    int p;
    switch (rbt) {
        case RBT_NONE: // 0
            return -1;  // this value is unused.
        case RBT_SYS: // 1
            return i;
        case RBT_CYC: // 2
            return (node_id + 1 + i) % n;
        default:
            p = rbt - RBT_SYSP;
            if (node_id <= p)
                return i; // node [0, p] help nodes [0, a) -- so node 0 has p preferred helpers.
            else
                return (i + 1) % n; // other nodes help [1, a+1) -- in particular, NOT node 0.
    }
}

/*
 * Returns the index [0, alpha) of the repair-chunk for node h_id helping f_id,
 * or -1 if it must be computed.
 */
int get_repair_chunk_idx(const int h_id, const int f_id,
        const int n, const int a, const RbtType rbt) {
    int idx;
    int p;
    switch (rbt) {
        case RBT_NONE: // 0
            return -1;
        case RBT_SYS: // 1
            if (f_id < a)
                return f_id;
            return -1;
        case RBT_CYC: // 2
            idx  = (f_id - h_id - 1 + n) % n;
            if (idx < a)
                return idx;
            return -1;
        default: // 3+ (SYSP)
            p = rbt - RBT_SYSP;
            if (h_id <= p && f_id < a)
                return f_id; // node [0, p] help nodes [0, a) -- so node 0 has p preferred helpers.
            else if (h_id > p && f_id >= 1 && f_id <= a)
                return f_id - 1; // other nodes help [1, a+1) -- in particular, NOT node 0.
            else
                return -1;
    }
}

/*
 * Helper to generate the RBT shuffling-matrix P
 * (extracting rows from node_mus).
 */
void get_rbt_P(int* P, const int node_id, const int* node_mus,
        const int n, const int a, const RbtType rbt) {
    for (int i = 0; i < a; i++) {
        int c = get_child_id(node_id, i, n, rbt);
        memcpy(P + i*a, node_mus + c*a, a*sizeof(int));
    }
}

// Generate P, P_inverse (+P_idx) matrices for all nodes.
void gen_rbt_data(int** all_P, int** all_P_inverses, int** nodes_P_idx,
        const int* node_mus, const int n, const int a, const int w,
        const RbtType rbt) {
    *nodes_P_idx = talloc(int, n);

    // TODO: if RBT_SYS, this can be done more efficiently
    // since all the P/P-inverses are the same.
    // (so use nodes_P_idx appropriately).

    *all_P = talloc(int, n*a*a);
    *all_P_inverses = talloc(int, n*a*a);
    int* P_scratch = talloc(int, a*a);
    for (int i = 0; i < n; i++) {
        // setup P for the helper i
        int* P = *all_P + i*(a*a);
        get_rbt_P(P, i, node_mus, n, a, rbt);

        // put inverses (all distinct) into P_inverses
        memcpy(P_scratch, P, a*a*sizeof(int));
        int* P_inv = *all_P_inverses + i*(a*a);
        invert_matrix(P_scratch, P_inv, a, w); // destroyes P_scratch

        (*nodes_P_idx)[i] = i;
    }
}

int* make_rbt_precoder_matrix(
        const int* all_P, const int* nodes_P_idx, const int* G,
        const int k, const int a, const int B, const int w) {
    // Compute the effective G for systematic nodes (after rbt postcoding).
    int* M = talloc(int, (k*a + B) * B); // the extra B*B is for augmentation.
    for (int i = 0; i < k; i++) {
        int* M_i = M + (a*B)*i;
        const int* G_i = G + (a*B)*i;
        const int* P_i = all_P + (a*a)*nodes_P_idx[i];
        matrix_multiply(M_i, P_i, G_i, a, a, a, B, w);
    }

    // Augment it as [M \\ I], then column-reduce.
    int* I_row = M + k*a*B;
    for (int i = 0; i < B; i++) {
      memset(I_row, 0, sizeof(int)*B);
      I_row[i] = 1;
      I_row += B;
    }
    cref(M, k*a + B, B, w);

    //jerasure_print_matrix(M, k*a + B, B, w);

    int* Q = talloc(int, B*B);
    int* Q_start = M + k*a*B;
    submat_cpy(
      Q_start, B,   // src
      B, B,         // B x B
      Q, B);        // dst
    free(M);

    return Q;
}

// P_inverses: concat of a x a rbt "unshuffling" matrices.
//             MUST be presistant... does not copy into plan.
// nodes_P_idx: k indices of the dc nodes, into P_inverses
GeneratorDAG* rbt_decode_dc(
        const int* P_inverses, const int* nodes_P_idx,
        const int k, const int a, const int w) {
    return rbt_shuffle(P_inverses, nodes_P_idx, k, a, w, "rbt_decode_root");
}

// Block-diagonal generator mtx for shuffle.
// (same operation as DAG version).
int* rbt_shuffle_mtx(
        const int* all_P, const int* nodes_P_idx,
        const int k, const int a, const int w) {
    int* G = (int*) calloc((k*a)*(k*a), sizeof(int));
    for (int i = 0; i < k; i++) {
        int* src = all_P + nodes_P_idx[i]*(a*a);
        int* dst = G + (i*a)*(k*a) + i*a;
        submat_cpy(src, a, a, a, dst, k*a);
    }
    return G;
}


// Returns a DAG that:
//  - splits a (ka)-length input into k alpha-length chunks
//  - multiplies vector i by P_i
GeneratorDAG* rbt_shuffle(
        const int* all_P, const int* nodes_P_idx,
        const int k, const int a, const int w,
        const char* stage_name) {
    GeneratorDAG* root = new_dag();
    root->stage_name = stage_name;
    root->stage_op = SHUFFLE;
    root->out_dim = k*a;
    root->I = range(0, k*a);
    root->num_deps = k;

    for (int i = 0; i < k; i++) {
        // first, extract the alpha symbols sent by dc-node i
        char* exName = (char*)malloc(40);
        sprintf(exName, "rbt extract (%d of %d)", i+1, k);
        GeneratorDAG* extractN = extract_from_dag(
                NULL, 1, k*a, // src
                0, i*a, // start (i, j)
                1, a, // num_rows, num_cols
                exName);

        // then apply P_inv[i] to these symbols
        GeneratorDAG* invertN = new_dag();
        invertN->stage_name = (char*)malloc(40);
        sprintf(invertN->stage_name, "rbt invert (%d of %d)", i+1, k);
        invertN->stage_op = TRANSFORM;
        invertN->out_dim = a;
        invertN->G = all_P + nodes_P_idx[i]*(a*a);
        invertN->num_deps = 1;
        invertN->deps[0] = extractN;

        root->deps[i] = invertN;
    }

    return root;
}


/*
 * Generates the mu_f helper-vector for an RBT-code,
 * given the original mu_f vector.
 */
int* rbt_mu_f(const int* orig_u_f, const int* P_inv,
        const int a, const int w) {
    int* mu_f = jerasure_matrix_multiply(orig_u_f, P_inv, 1, a, a, a, w);
    return mu_f;
}
