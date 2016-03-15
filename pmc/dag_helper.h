#ifndef DAG_HELPER_H
#define DAG_HELPER_H

#include "product_matrix.h"
#include <stdint.h>

#define TAG_DONE 0
#define TAG_TODO ((void*)(-1))

// Time required to process top-level DAG (excluding dependencies)
double get_exclusive_time(GeneratorDAG* dag);

void print_dag(GeneratorDAG* dag);

// Returns a new, empty dag.
// ALWAYS use this to alloc a DAG, because it
// initializes the DONE and TAG fields to zero.
GeneratorDAG* new_dag();

GeneratorDAG* concat_dag(const GeneratorDAG* first, GeneratorDAG* second);

int* generator_mtx_from_dag(GeneratorDAG* dag, int in_dim, int w);

// Merges root and root->deps[0] into a single transform.
// MUTATES root->deps[0]
GeneratorDAG* merge_top2_transforms(GeneratorDAG* root, int w);

void reset_done(GeneratorDAG* dag);

void reset_tags(GeneratorDAG* dag, void* value, const int clear_done);

int* range(int start_inc, int end_ex);

// helper-dag: esentially a no-op, useful placeholder for composing
// the message matrix from which to extract msg vector.
GeneratorDAG* noop_dag(int out_dim);

// Wraps a simple generator-matrix multiplication as a dag stage.
// handles SRC=NULL appropriately
GeneratorDAG* make_transform_dag(const GeneratorDAG* src, const int* G,
        const int out_dim, const char* name);

// Unwraps a decoded message matrix
// (stored in M_dag with structure specified in M),
// into a vector of the unique symbols: [m0, m1, m2,...]
GeneratorDAG* msg_unwrap_dag(GeneratorDAG* M_dag, int* M, int d, int a, int B);

// returns the SHUFFLE-node for extracting a submatrix.
// handles SRC=NULL properly (--> num_deps = 0)
GeneratorDAG* extract_from_dag(GeneratorDAG* src, int src_rows, int src_cols, int i, int j, int nr, int nc, char* name);

// A*src, where A: n x k, src: k x m
GeneratorDAG* left_mult_dag(int* A, GeneratorDAG* src, int n, int k, int m, char* name);

#define sub_dag add_dag
GeneratorDAG* add_dag(GeneratorDAG* A, GeneratorDAG* B, int r, int c, char* name);

// returns the dag SHUFFLE-node that hstacks [A | B]
GeneratorDAG* hstack_dag(GeneratorDAG* A, GeneratorDAG* B, int A_cols, int B_cols, int rows, char* name);

// returns the dag SHUFFLE-node that vstacks [A \\ B]
GeneratorDAG* vstack_dag(GeneratorDAG* A, GeneratorDAG* B, int A_rows, int B_rows, int cols, char* name);

GeneratorDAG* transpose_dag(GeneratorDAG* src, int src_rows, int src_cols);

// converts a plan into a DAG, using input_dag as the base input
GeneratorDAG* dag_from_plan(GeneratorPlan* plan, GeneratorDAG* input_dag);

#endif
