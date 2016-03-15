#include "product_matrix.h"
#include <string.h>
#include "matrix.h"
#include "MBR_MSR.h"
#include "../jerasure-kmg/include/jerasure.h"
#include "math.h"
#include <sys/time.h>
#include <time.h>

#include "dag_helper.h"

#define talloc(type, num) (type *) malloc(sizeof(type)*(num))

void make_systematic(int* G, int n, int a, int B, int w)
{
  cref(G, (n*a), B, w);
}

/*
  Returns the (non-systematic) generator mtx for PMC specified by Psi and M
*/
int* generator_for_PMC(int* psi, int* M, int n, int d, int a, int B)
{
  //M is (d x a)
  int* G = talloc(int, (n*a)*B); // G is (na) x B
  memset(G, 0, sizeof(int)*(n*a)*B);

  int i, j, s;
  for(i = 0; i < n; i++)
    for (j = 0; j < a; j++) { //the jth stored symbol of the ith node
      //row (i*a + j) of G
      int r = i*a + j;

      for (s = 0; s < d; s++) {
        int e = M[s*a + j]; // M[s, j]
        if(e > 0)
          G[r*B + (e-1)] = psi[i*d + s]; //G[r, (e-1) = psi[i, s]
      }
    }

  return G;
}


/*
 * Returns the matrix G such that G * vec(X) = vec(AX)
 * where vec(.) is the matrix-vectorization operator.
 *
 * A : n x k
 * X : k x m
 */
int* left_mult_G(int* A, int n, int k, int m)
{
  int* G = talloc(int, (n*m)*(k*m)); // G is (nm) x (km)
  memset(G, 0, sizeof(int)*(n*m)*(k*m));

  int i, j, s;
  for(i = 0; i < n; i++)
    for (j = 0; j < m; j++) { // (i,j)th element in matrix product AX
      //row (i*m + j) of G
      int r = i*m + j;

      for (s = 0; s < k; s++) {
          G[r*(k*m) + (s*m + j)] = A[i*k + s]; //G[r, sm + j] = A[i, s]
      }
    }

  return G;
}

void right_mult_G_into(int* G, int* A, int n, int k, int m)
{
  memset(G, 0, sizeof(int)*(n*m)*(n*k));
  int i, j, s;
  for(i = 0; i < n; i++)
    for (j = 0; j < m; j++) { // (i,j)th element in matrix product XA
      //row (i*m + j) of G
      int r = i*m + j;

      for (s = 0; s < k; s++) {
          G[r*(n*k) + (i*k + s)] = A[s*m + j]; //G[r, ik + s] = A[s, j]
      }
    }
}

/*
 * Returns the matrix G such that G * vec(X) = vec(XA)
 * where vec(.) is the matrix-vectorization operator.
 *
 * X : n x k
 * A : k x m
 */
int* right_mult_G(int* A, int n, int k, int m) {
    int* G = talloc(int, (n*m)*(n*k)); // G is (nm) x (nk)
    right_mult_G_into(G, A, n, k, m);
    return G;
}

int* apply_plan_single(GeneratorPlan* plan, int* msg, int w) {
    int* input = msg;
    int* output;

    int** D = plan->D;
    int* dim = plan->dim;

    //printf("starting plan...\n");
    while(*D) {
        int inDim = *dim;
        int outDim = *(++dim);

        //printf("applying D (%d x %d):\n", outDim, inDim);
        //jerasure_print_matrix(*D, outDim, inDim, w);
        output = jerasure_matrix_multiply(*D, input, outDim, inDim, inDim, 1, w);

        //printf("output:\n");
        //jerasure_print_matrix(output, 1, outDim, w);

        if(input != msg)
            free(input);
        input = output;
        D++;
    }
    return output;
}

void apply_dag_striped_buff(GeneratorDAG* dag, char** base_input_ptrs, char** output_ptrs,
        int base_input_dim, int chunk_size, int w, int buffsize) {
    dag->start_t = clock();
    char** input_ptrs;

    if (dag->num_deps == 0) {
        input_ptrs = base_input_ptrs;
        dag->in_dim = base_input_dim;
    } else {
        // generate dependencies
        dag->in_dim = 0;
        for (int i = 0; i < dag->num_deps; i++) {
            GeneratorDAG* dep = dag->deps[i];
            dag->in_dim += dep->out_dim;
            if (! dep->done) {
                char** dep_output_ptrs = talloc(char*, dep->out_dim);

                apply_dag_striped(dep, base_input_ptrs, dep_output_ptrs,
                        base_input_dim, chunk_size, w);
                dep->tag = dep_output_ptrs; // store the output ptrs in the tag
            }
        }

        // setup input_ptrs sequentially
        int n = 0;
        input_ptrs = talloc(char*, dag->in_dim);
        for (int i = 0; i < dag->num_deps; i++)
            for (int j = 0; j < dag->deps[i]->out_dim; j++)
                input_ptrs[n++] = ((char**)dag->deps[i]->tag)[j];
    }

    if (dag->stage_op == SHUFFLE) {
        for (int i = 0; i < dag->out_dim; i++)
            output_ptrs[i] = input_ptrs[dag->I[i]];
    } else {
        for (int i = 0; i < dag->out_dim; i++) {
            aligned_malloc_pmc(output_ptrs + i, chunk_size);
        }
        matrix_data_product_buff(dag->G, dag->in_dim, dag->out_dim,
                input_ptrs, output_ptrs, chunk_size,
                w, buffsize);
    }
    dag->done = 1;
    dag->end_t = clock();
}

/*
 * Executes DAG with base_input_ptrs.
 * Each (of base_input_dim) src_ptr[i] points to chunk_size bytes of data.
 * Sets output_ptrs to point to result (allocating space).
 */
void apply_dag_striped(GeneratorDAG* dag, char** base_input_ptrs, char** output_ptrs,
        int base_input_dim, int chunk_size, int w) {
    dag->start_t = clock();
    char** input_ptrs;

    if (dag->num_deps == 0) {
        input_ptrs = base_input_ptrs;
        dag->in_dim = base_input_dim;
    } else {
        // generate dependencies
        dag->in_dim = 0;
        for (int i = 0; i < dag->num_deps; i++) {
            GeneratorDAG* dep = dag->deps[i];
            dag->in_dim += dep->out_dim;
            if (! dep->done) {
                char** dep_output_ptrs = talloc(char*, dep->out_dim);

                apply_dag_striped(dep, base_input_ptrs, dep_output_ptrs,
                        base_input_dim, chunk_size, w);
                dep->tag = dep_output_ptrs; // store the output ptrs in the tag
            }
        }

        // setup input_ptrs sequentially
        int n = 0;
        input_ptrs = talloc(char*, dag->in_dim);
        for (int i = 0; i < dag->num_deps; i++)
            for (int j = 0; j < dag->deps[i]->out_dim; j++)
                input_ptrs[n++] = ((char**)dag->deps[i]->tag)[j];
    }

    if (dag->stage_op == SHUFFLE) {
        for (int i = 0; i < dag->out_dim; i++)
            output_ptrs[i] = input_ptrs[dag->I[i]];
    } else {
        int* G_i = dag->G;
        for (int i = 0; i < dag->out_dim; i++) {
            //output_ptrs[i] = (char*) malloc(chunk_size);
            aligned_malloc_pmc(output_ptrs + i, chunk_size);
            row_data_product(G_i, input_ptrs, output_ptrs[i], dag->in_dim, chunk_size, w);
            G_i += dag->in_dim;
        }
    }
    dag->done = 1;
    dag->end_t = clock();
    // debugging
#ifdef DEBUG
    printf("STAGE DONE: %s\n", dag->stage_name);

    /*if (strcmp(dag->stage_name, "M") == 0) {*/
        /*for (int i = 0; i < dag->out_dim; i++) {*/
            /*printf("chunk %d:", i);*/
            /*for (int b = 0; b < chunk_size; b++)*/
                /*putchar(output_ptrs[i][b]);*/
            /*printf("\n");*/
        /*}*/
    /*}*/
#endif
}

/*
 * Applies plan across src_ptrs.
 * Each src_ptr[i] points to chunk_size bytes of data.
 */
void apply_plan_striped(GeneratorPlan* plan, char** src_ptrs, char** dst_ptrs, int chunk_size, int w) {
    char** input_ptrs = src_ptrs;
    char** output_ptrs;

    int** D = plan->D;
    int* dim = plan->dim;

    //printf("starting plan...\n");
    while(*D) {
        int inDim = *dim;
        int outDim = *(++dim);

        if (*(D+1)) {
            output_ptrs = talloc(char*, outDim);
            for (int i = 0; i < outDim; i++)
                output_ptrs[i] = (char*) malloc(chunk_size);
        } else {
            output_ptrs = dst_ptrs;
        }

        //printf("applying D (%d x %d):\n", outDim, inDim);
        int* D_i = *D;
        for (int i = 0; i < outDim; i++) {
            row_data_product(D_i, input_ptrs, output_ptrs[i], inDim, chunk_size, w);
            D_i += inDim;
        }

        if(input_ptrs != src_ptrs) {
            for (int i = 0; i < inDim; i++)
                free(input_ptrs[i]);
            free(input_ptrs);
        }

        input_ptrs = output_ptrs;
        D++;
    }
}

void apply_data_collect_plan_striped(GeneratorPlan* plan, char** chunk_ptrs, char* dst, int B, int msg_size, int w) {
    int chunk_size = msg_size / B;

    char** dst_ptrs = talloc(char*, B);
    for (int i = 0; i < B; i++) {
        dst_ptrs[i] = dst + i*chunk_size;
    }

    apply_plan_striped(plan, chunk_ptrs, dst_ptrs, chunk_size, w);
}

void apply_data_collect_dag_striped(GeneratorDAG* dag, char** chunk_ptrs, char* dst, int k, int a, int B, int msg_size, int w) {
    //TODO: avoid so much data copy...
    int chunk_size = msg_size / B;

    char** dst_ptrs = talloc(char*, B);
    apply_dag_striped(dag, chunk_ptrs, dst_ptrs, k*a, chunk_size, w);

    for (int i = 0; i < B; i++) {
        memcpy(dst + i*chunk_size, dst_ptrs[i], chunk_size);
    }
}

/*
 * Convenience function to set up B pointers into MSG
 * (splitting the message into (msg_size / B)-sized chunks for each symbol)
 */
char** init_data_ptrs(const char* msg, const int B, const int msg_size) {
  //int sym_size = w / 8;
  //int stripe_size = B * sym_size;
  //int num_stripes = msg_size / stripe_size;

  if (msg_size % B != 0) {
    fprintf(stderr, "ERROR: msg_size is not a multiple of B\n");
    exit(1);
  }
  if ((msg_size / B) % sizeof(long) != 0) {
    fprintf(stderr, "ERROR: (msg_size / B) = (%d / %d) is not a multiple of sizeof(long)\n", msg_size, B);
    exit(1);
  }

  /* Setup B data pointers into msg. */
  char** data_ptrs = talloc(char*, B);
  int data_size = msg_size / B;
  for (int i = 0; i < B; i++) {
    data_ptrs[i] = msg + i*data_size; // <==> ... + i*(num_stripes * sym_size)
  }
  return data_ptrs;
}

/*
  Encodes for a specific node i (generator Gi).
  Splits msg bytes into B-symbol-sized stripes,
  takes the product Gi * x for each stripe x,
  and puts the result of the (j-th row of Gi) * (all stripes x) into dst_ptrs[j]
  That is, dst_ptrs[j] contains the j-th symbol of all stripes. (sizes: dst_ptrs[alpha][num_stripes*sym_size])

  msg_size in bytes.
*/
void encode_striped(char* msg, int* Gi, char** dst_ptrs, int a, int B, int msg_size, int w)
{
  char** data_ptrs = init_data_ptrs(msg, B, msg_size);
  int data_size = msg_size / B;

  /* Encode all stripes for every one of a=alpha symbols. */
  int* Gi_j = Gi; //j-th row of Gi --> j-th symbol.
  for (int j = 0; j < a; j++) {
    row_data_product(Gi_j, data_ptrs, dst_ptrs[j], B, data_size, w);
    Gi_j += B;
  }
}

/*
  Takes in pointers to (k*a) chunks of recieved data (each chunk msg_size / B bytes),
  and writes the recovered msg_size bytes of message to dst.
*/
void data_collect_striped(char** chunk_ptrs, int* D, char* dst, int k, int a, int B, int msg_size, int w)
{
  int chunk_size = msg_size / B;

  int* D_i = D;
  for (int i = 0; i < B; i++) {
    /* Decode the i-th symbol of all stripes. */
    row_data_product(D_i, chunk_ptrs, dst, (k*a), chunk_size, w);
    D_i += (k*a);
    dst += chunk_size;
  }
}

/*
  Takes in pointers to (d) chunks of recieved data (each chunk msg_size / B bytes),
  and writes the i-th repaired chunk to dst_ptrs[i]. (of alpha total repaired chunks).
*/
void repair_striped(char** chunk_ptrs, int* A, char** dst_ptrs, int d, int a, int B, int msg_size, int w)
{
  int chunk_size = msg_size / B;

  int* A_i = A;
  for (int i = 0; i < a; i++) {
    row_data_product(A_i, chunk_ptrs, dst_ptrs[i], d, chunk_size, w);
    A_i += d;
  }
}

void helper_symbol_striped(char** chunk_ptrs, int* u_f, char* dst, int a, int B, int msg_size, int w) {
  int chunk_size = msg_size / B;
  row_data_product(u_f, chunk_ptrs, dst, a, chunk_size, w);
}

int* encode(int* msg, int* G, int n, int a, int B, int w)
{
  return jerasure_matrix_multiply(G, msg, n*a, B, B, 1, w);
}

int* helper_symbol(int* u_f, int* y, int a, int w)
{
  return jerasure_matrix_multiply(u_f, y, 1, a, a, 1, w);
}

int* data_collect(int* D, int* yDC, int k, int a, int B, int w)
{
  return jerasure_matrix_multiply(D, yDC, B, k*a, k*a, 1, w);
}

int* repair(int* A, int* yR, int d, int a, int w)
{
  return jerasure_matrix_multiply(A, yR, a, d, d, 1, w);
}


/*
  Returns the (B x ka) data-collection matrix D
  such that if y (ka x 1) is the recieved vector from k nodes,
  then Dy (B x 1) is the recovered message.
*/
int* make_data_collect_matrix(int* G, int k, int a, int B, int* nodes, int w)
{
  int gi_r = k*a;
  int gi_c = B + k*a;
  int* GI = talloc(int, gi_r * gi_c); // [[G]_k | I]

  /* Copy the appropriate generator rows from participating nodes. */
  for (int i = 0; i < k; i++) {
    int* G_ni = G + (nodes[i]*a) * B;
    int* GI_targ = GI + (i*a) * gi_c;

    submat_cpy(
      G_ni, B,        // src
      a, B,           // a x B
      GI_targ, gi_c); // dst
  }

  /* Set up identity (augmented on right). */
  int* I_row = GI + B;
  for (int i = 0; i < k*a; i++) {
    memset(I_row, 0, sizeof(int)*(k*a));
    I_row[i] = 1;
    I_row += gi_c;
  }


  /* Row-reduce, and put the top B rows of the augmented portion in D. */

  rref_aug(GI, gi_r, gi_c, B, w);
  //jerasure_print_matrix(GI, k*a, B + k*a, w);

  int* D = talloc(int, B * (k*a));

  int* D_start = GI + B;
  submat_cpy(
    D_start, gi_c,  // src
    B, k*a,         // B x (ka)
    D, k*a);        // dst
  free(GI);

  return D;
}

int* make_repair_matrix(int* G, int d, int a, int B, int f, int* u_f, int* helper_nodes, int w)
{
  int frows = B;
  int fcols = (d + a);
  int* F = talloc(int, frows * fcols);

  /* Compute W.  */
  int* W = talloc(int, d * B);  //w_i = talloc(int, B);
  for (int i = 0; i < d; i++) {
    int* G_i = G + (helper_nodes[i]*a) * B;
    int* w_i = W + i*B;
    matrix_multiply(w_i, u_f, G_i, 1, a, a, B, w); // into w_i
  }

  transpose_block(W, F, d, B, B, fcols); // put W^T into F
  free(W);

  int* G_f = G + (f*a) * B;
  transpose_block(G_f, F + d, a, B, B, fcols); // put G_f^T into F (augmenting W^T)

  rref_aug(F, frows, fcols, d, w);
  //jerasure_print_matrix(F, frows, fcols, w);

  int* A = talloc(int, a * d);
  transpose_block(F + d, A, d, a, fcols, d); // pull out first d rows of augmented part as A^T.
  free(F);

  return A;
}
