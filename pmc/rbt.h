#ifndef RBT_H
#define RBT_H

#define talloc(type, num) (type *) malloc(sizeof(type)*(num))

#include "product_matrix.h"
#include "matrix.h"

typedef enum {RBT_NONE = 0, RBT_SYS, RBT_CYC, RBT_SYSP} RbtType;

int get_child_id(const int node_id, const int i, const int n,
        const RbtType rbt);

int get_repair_chunk_idx(const int h_id, const int f_id,
        const int n, const int a, const RbtType rbt);

/*
 * Helper to generate the RBT shuffling-matrix P
 * (extracting rows from node_mus).
 */
void get_rbt_P(int* P, const int node_id, const int* node_mus,
        const int n, const int a, const RbtType rbt);

void gen_rbt_data(int** P, int** P_inverses, int** nodes_P_idx,
        const int* node_mus, const int n, const int a, const int w,
        const RbtType rbt);

/* For RBT Encoding. */

int* make_rbt_precoder_matrix(
        const int* all_P, const int* nodes_P_idx, const int* G,
        const int k, const int a, const int B, const int w);

/* For RBT Data-Collecting. */

// Block-diagonal generator mtx for shuffle.
// (same operation as DAG version).
int* rbt_shuffle_mtx(
        const int* all_P, const int* nodes_P_idx,
        const int k, const int a, const int w);

GeneratorDAG* rbt_shuffle(
        const int* all_P, const int* nodes_P_idx,
        const int k, const int a, const int w,
        const char* stage_name);

GeneratorDAG* rbt_decode_dc(
        const int* P_inverses, const int* nodes_P_idx,
        const int k, const int a, const int w);


/* For RBT Repair. */

int* rbt_mu_f(const int* orig_u_f, const int* P_inv,
        const int a, const int w);

int get_repair_chunk_idx(const int h_id, const int f_id, const int n, const int a, const RbtType rbt);

#endif
