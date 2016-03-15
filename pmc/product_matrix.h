#ifndef PRODUCT_MATRIX_H
#define PRODUCT_MATRIX_H

#include <sys/time.h>
#include <time.h>

#include <stdlib.h>
#define aligned_malloc_pmc(handle, bytes) posix_memalign(handle, 16, bytes);

typedef struct {
    int* D[50];
    int dim[51];
} GeneratorPlan;

typedef enum {SHUFFLE, TRANSFORM} StageOp;
typedef struct GeneratorDAG {
    StageOp stage_op;
    char* stage_name;
    int done;
    clock_t start_t, end_t;
    void* tag;

    int out_dim;
    int in_dim;

    int num_deps; // num dependencies. if 0, then operates on message
    struct GeneratorDAG* deps[50];

    int* G; // Generator matrix: only for transform
    int* I; // Shuffle index-representation: only for shuffle (not strictly neccesary to keep seperate)
    // output[i] = input[I[i]]

} GeneratorDAG;

/*
  Returns the (non-systematic) generator mtx for PMC specified by Psi and M
*/
int* generator_for_PMC(int* psi, int* M, int n, int d, int a, int B);

/*
 * Returns the matrix G such that G * vec(X) = vec(AX)
 * where vec(.) is the matrix-vectorization operator.
 *
 * A : n x k
 * X : k x m
 */
int* left_mult_G(int* A, int n, int k, int m);

/*
 * Returns the matrix G such that G * vec(X) = vec(XA)
 * where vec(.) is the matrix-vectorization operator.
 *
 * X : n x k
 * A : k x m
 */
int* right_mult_G(int* A, int n, int k, int m);

void right_mult_G_into(int* G, int* A, int n, int k, int m);

void make_systematic(int* G, int n, int a, int B, int w);

int* apply_plan_single(GeneratorPlan* plan, int* msg, int w);

void apply_dag_striped_buff(GeneratorDAG* dag, char** base_input_ptrs, char** output_ptrs,
        int base_input_dim, int chunk_size, int w, int buffsize);

void apply_dag_striped(GeneratorDAG* dag, char** base_input_ptrs, char** output_ptrs, int base_input_dim, int chunk_size, int w);

/*
 * Applies plan across src_ptrs.
 * Each src_ptr[i] points to chunk_size bytes of data.
 */
void apply_plan_striped(GeneratorPlan* plan, char** src_ptrs, char** dst_ptrs, int chunk_size, int w);

/*
 * Applies a data-collection plan to chunk_ptrs (length B), and writes the
 * decoded message to dst.
 * Note: msg_size = B * chunk_size
 */
void apply_data_collect_plan_striped(GeneratorPlan* plan, char** chunk_ptrs, char* dst, int B, int msg_size, int w);


void apply_data_collect_dag_striped(GeneratorDAG* dag, char** chunk_ptrs, char* dst, int k, int a, int B, int msg_size, int w);

char** init_data_ptrs(const char* msg, const int B, const int msg_size);

/*
  Encodes for a specific node i (generator Gi).
  Splits msg bytes into B-symbol-sized stripes,
  takes the product Gi * x for each stripe x,
  and puts the result of the (j-th row of Gi) * (all stripes x) into dst_ptrs[j]
  That is, dst_ptrs[j] contains the j-th symbol of all stripes. (sizes: dst_ptrs[alpha][num_stripes*sym_size])

  msg_size in bytes.
*/
void encode_striped(char* msg, int* Gi, char** dst_ptrs, int a, int B, int msg_size, int w);

/*
  Takes in pointers to (k*a) chunks of recieved data (each chunk msg_size / B bytes),
  and writes the recovered msg_size bytes of message to dst.
*/
void data_collect_striped(char** chunk_ptrs, int* D, char* dst, int k, int a, int B, int msg_size, int w);

/*
  Takes in pointers to (d) chunks of recieved data (each chunk msg_size / B bytes),
  and writes the i-th repaired chunk to dst_ptrs[i]. (of alpha total repaired chunks).
*/
void repair_striped(char** chunk_ptrs, int* A, char** dst_ptrs, int d, int a, int B, int msg_size, int w);

void helper_symbol_striped(char** chunk_ptrs, int* u_f, char* dst, int a, int B, int msg_size, int w);

int* encode(int* msg, int* G, int n, int a, int B, int w);

int* helper_symbol(int* u_f, int* y, int a, int w);

int* data_collect(int* D, int* yDC, int k, int a, int B, int w);

int* repair(int* A, int* yR, int d, int a, int w);

/*
  Returns the (B x ka) data-collection matrix D
  such that if y (ka x 1) is the recieved vector from k nodes,
  then Dy (B x 1) is the recovered message.
*/
int* make_data_collect_matrix(int* G, int k, int a, int B, int* nodes, int w);

int* make_repair_matrix(int* G, int d, int a, int B, int f, int* u_f, int* helper_nodes, int w);

#endif
