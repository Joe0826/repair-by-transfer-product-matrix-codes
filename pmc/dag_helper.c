#include "dag_helper.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>

#include "matrix.h"

#define talloc(type, num) (type *) malloc(sizeof(type)*(num))

// Time required to process top-level DAG (excluding dependencies)
double get_exclusive_time(GeneratorDAG* dag) {
    clock_t deps_done_t = dag->start_t;
    for (int i = 0; i < dag->num_deps; i++)
        deps_done_t = fmax(deps_done_t, dag->deps[i]->end_t);
    return ((double)(dag->end_t - deps_done_t)/CLOCKS_PER_SEC);
}

static void print_dag_lev(GeneratorDAG* dag, int level) {
    for (int t = 0; t < level; t++)
        printf("|  ");
    double dur = ((double)(dag->end_t - dag->start_t) / CLOCKS_PER_SEC);
    double exdur = get_exclusive_time(dag);
    printf("%s [total: %.5f, excl: %.5f]\n", dag->stage_name, dur, exdur);

    for (int i = 0; i < dag->num_deps; i++)
        print_dag_lev(dag->deps[i], level+1);
}
void print_dag(GeneratorDAG* dag) {
    print_dag_lev(dag, 0);
}

GeneratorDAG* new_dag() {
    return (GeneratorDAG*) calloc(1, sizeof(GeneratorDAG));
}

static GeneratorDAG* concat_dag_helper(const GeneratorDAG* first,
    GeneratorDAG* second) {
    GeneratorDAG* root = second;

    if (root->tag != TAG_DONE) {
        if (root->num_deps == 0) { // would act on base input...
            root->num_deps = 1;
            root->deps[0] = first; // point to FIRST instead
        } else {
            for (int i = 0; i < root->num_deps; i++) {
                GeneratorDAG* dep = root->deps[i];
                concat_dag_helper(first, dep);
            }
        }
        root->tag = TAG_DONE;
    }
    return root;
}

// returns the DAG that : applies FIRST, then applies SECOND
// ie, the depnedency dag is: SECOND --> FIRST --> <base>
// MUTATES second into the final dag.
// (also messes up the tags)
// TODO(preetum,jingyan): add tail-pointers to DAGs, so this can be
// done without traversal.
GeneratorDAG* concat_dag(const GeneratorDAG* first, GeneratorDAG* second) {
    if (first == NULL)
        return second;
    else if (second == NULL)
        return first;

    reset_tags(second, TAG_TODO, 0);
    return concat_dag_helper(first, second);
}

static void compute_in_dim(GeneratorDAG* dag) {
    dag->in_dim = 0;
    for (int i = 0; i < dag->num_deps; i++) {
        GeneratorDAG* dep = dag->deps[i];
        dag->in_dim += dep->out_dim;
    }
}

int* generator_mtx_from_dag(GeneratorDAG* dag, int in_dim, int w) {
    if (w != 8) {
        fprintf(stderr, "ERROR: only w=8 supported.\n");
        exit(1);
    }
    int out_dim = dag->out_dim;
    int* G = (int*)calloc(out_dim*in_dim, sizeof(int));
    char** e_i = (char**)malloc(sizeof(char*)*in_dim);
    for (int i = 0; i < in_dim; i++)
        aligned_malloc_pmc(e_i + i, 4);

    char** output_ptrs = (char**)malloc(sizeof(char*)*out_dim);
    for (int i = 0; i < in_dim; i++) {
        // e_i = i-th basis vector
        for (int j = 0; j < in_dim; j++)
            e_i[j][0] = (i==j)? 1 : 0;

        apply_dag_striped_buff(dag, e_i, output_ptrs, in_dim, 1, w, 1000);
        reset_done(dag);

        for (int k = 0; k < out_dim; k++) {
            unsigned char val = output_ptrs[k][0];
            G[i + k*in_dim] = val;
        }
    }
    return G;
}

// Merges root and root->deps[0] into a single transform.
// MUTATES root->deps[0]
GeneratorDAG* merge_top2_transforms(GeneratorDAG* root, int w) {
    GeneratorDAG* dep = root->deps[0];
    if (root->stage_op != TRANSFORM || root->num_deps != 1
            || dep->stage_op != TRANSFORM) {
        fprintf(stderr, "ERROR: merge dags malformed\n");
        return NULL;
    }
    int dep_out_dim = dep->out_dim;
    char* mname = (char*)malloc(100);
    sprintf(mname, "%s <merged with> %s", root->stage_name, dep->stage_name);
    GeneratorDAG* merged = dep;
    merged->stage_name = mname;
    merged->out_dim = root->out_dim;
    compute_in_dim(merged);

    int* G = (int*) malloc(sizeof(int)*(merged->out_dim)*(merged->in_dim));
    matrix_multiply(G, root->G, dep->G, merged->out_dim, dep_out_dim,
            dep_out_dim, merged->in_dim, w);

    merged->G = G;
    return merged;
}

void reset_done(GeneratorDAG* dag) {
    dag->done = 0;
    for (int i = 0; i < dag->num_deps; i++) {
        GeneratorDAG* dep = dag->deps[i];
        if (dep->done)
            reset_done(dep);
    }
}

// resets all tags in dag to value.
// FAILS if any tags were already set to value... TODO(preetum): fix this
void reset_tags(GeneratorDAG* dag, void* value, const int clear_done) {
    if (clear_done)
        dag->done = 0;
    if (dag->tag != (void*) value) {
        for (int i = 0; i < dag->num_deps; i++) {
            GeneratorDAG* dep = dag->deps[i];
            reset_tags(dep, value, clear_done);
        }
        dag->tag = (void*) value;
    }
}

// returns an array [start_inc, end_ex)
int* range(int start_inc, int end_ex) {
    int* I = talloc(int, end_ex - start_inc);
    for (int i = start_inc; i < end_ex; i++)
        I[i - start_inc] = i;
    return I;
}

// helper-dag: esentially a no-op, useful placeholder for composing
// the message matrix from which to extract msg vector.
GeneratorDAG* noop_dag(int out_dim) {
    GeneratorDAG* NOOP_DAG = new_dag();
    NOOP_DAG->stage_name = "noop";
    NOOP_DAG->stage_op = SHUFFLE;
    NOOP_DAG->num_deps = 0;
    NOOP_DAG->out_dim = out_dim;
    NOOP_DAG->I = (int*)calloc(out_dim, sizeof(int));
    return NOOP_DAG;
}

// Wraps a simple generator-matrix multiplication as a dag stage.
// handles SRC=NULL appropriately
GeneratorDAG* make_transform_dag(const GeneratorDAG* src, const int* G,
        const int out_dim, const char* name) {
    GeneratorDAG* dag = new_dag();
    dag->stage_name = name;
    dag->stage_op = TRANSFORM;
    dag->out_dim = out_dim;
    if (src) {
        dag->deps[0] = src;
        dag->num_deps = 1;
    } else {
        dag->num_deps = 0;
    }
    dag->G = G;
    return dag;
}

GeneratorDAG* msg_unwrap_dag(GeneratorDAG* M_dag, int* M, int d, int a, int B) {
    int* I = talloc(int, B);
    for (int i = 0; i < d*a; i++) {
        int e = M[i] - 1;
        if (e >= 0)
            I[e] = i;
    }

    GeneratorDAG* msg = new_dag();
    msg->stage_name = "final unwrap";
    msg->stage_op = SHUFFLE;
    msg->out_dim = B;
    msg->I = I;
    msg->deps[0] = M_dag;
    msg->num_deps = 1;
    return msg;
}

// returns the SHUFFLE-node for extracting a submatrix.
// handles SRC=NULL properly (--> num_deps = 0)
GeneratorDAG* extract_from_dag(GeneratorDAG* src, int src_rows, int src_cols, int i, int j, int nr, int nc, char* name) {

    int* I = (int*) malloc(nr*nc*sizeof(int));
    for (int r = 0; r < nr; r++)
        for (int c = 0; c < nc; c++)
            I[r*nc + c] = (i+r)*src_cols + (j+c);

    GeneratorDAG* D = new_dag();
    D->stage_name = name;
    D->stage_op = SHUFFLE;
    D->out_dim = nr*nc;
    D->I = I;
    if (src) {
        D->deps[0] = src;
        D->num_deps = 1;
    } else {
        D->num_deps = 0;
    }
    return D;
}

// A*src, where A: n x k, src: k x m
GeneratorDAG* left_mult_dag(int* A, GeneratorDAG* src, int n, int k, int m, char* name) {
    GeneratorDAG* D = new_dag();
    D->stage_name = name;
    D->stage_op = TRANSFORM;
    D->out_dim = n*m;
    D->G = left_mult_G(A, n, k, m);
    D->deps[0] = src;
    D->num_deps = 1;
    return D;
}

GeneratorDAG* add_dag(GeneratorDAG* A, GeneratorDAG* B, int r, int c, char* name) {


    int* G = (int*) calloc((r*c)*(2*r*c), sizeof(int));
    for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++)
            G[(i*c + j)*(2*r*c) + i*c + j] = G[(i*c + j)*(2*r*c) + i*c + j + r*c] = 1;

    GeneratorDAG* D = new_dag();
    D->stage_name = name;
    D->stage_op = TRANSFORM;
    D->out_dim = r*c;
    D->G = G;
    D->deps[0] = A;
    D->deps[1] = B;
    D->num_deps = 2;
    return D;
}

// returns the dag SHUFFLE-node that hstacks [A | B]
GeneratorDAG* hstack_dag(GeneratorDAG* A, GeneratorDAG* B, int A_cols, int B_cols, int rows, char* name) {

    int* I = talloc(int, rows*(A_cols + B_cols));
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < A_cols + B_cols; j++) {
            int* idx = I + i*(A_cols + B_cols) + j;
            if (j < A_cols)
                *idx = i*A_cols + j;
            else
                *idx = rows*A_cols + i*B_cols + (j-A_cols);
        }

    GeneratorDAG* D = new_dag();
    D->stage_name = name;
    D->stage_op = SHUFFLE;
    D->out_dim = rows * (A_cols + B_cols);
    D->I = I;
    D->deps[0] = A;
    D->deps[1] = B;
    D->num_deps = 2;
    return D;
}

// returns the dag SHUFFLE-node that vstacks [A \\ B]
GeneratorDAG* vstack_dag(GeneratorDAG* A, GeneratorDAG* B, int A_rows, int B_rows, int cols, char* name) {
    int* I = talloc(int, (A_rows + B_rows)*cols);
    for (int i = 0; i < (A_rows + B_rows)*cols; i++)
            I[i] = i;

    GeneratorDAG* D = new_dag();
    D->stage_name = name;
    D->stage_op = SHUFFLE;
    D->out_dim = (A_rows + B_rows)*cols;
    D->I = I;
    D->deps[0] = A;
    D->deps[1] = B;
    D->num_deps = 2;
    return D;
}

GeneratorDAG* transpose_dag(GeneratorDAG* src, int src_rows, int src_cols) {
    int* I = talloc(int, src_cols*src_rows);
    for (int i = 0; i < src_cols; i++)
        for (int j = 0; j < src_rows; j++)
            I[i*src_rows + j] = j*src_cols + i;

    GeneratorDAG* D = new_dag();
    D->stage_name = "transpose";
    D->stage_op = SHUFFLE;
    D->out_dim = src_rows*src_cols;
    D->I = I;
    D->deps[0] = src;
    D->num_deps = 1;
    return D;
}

// converts a plan into a DAG, using input_dag as the base input
GeneratorDAG* dag_from_plan(GeneratorPlan* plan, GeneratorDAG* input_dag) {
    GeneratorDAG* dag;

    int** D = plan->D;
    int* dim = plan->dim;
    int i = 0;
    while(*D) {
        int inDim = *dim;
        int outDim = *(++dim);

        dag = new_dag();
        dag->stage_name = (char*)malloc(20);
        sprintf(dag->stage_name, "plan (stage %d)", i++);
        dag->stage_op = TRANSFORM;
        dag->out_dim = outDim;
        if (input_dag == NULL) {
            dag->num_deps = 0;
        } else {
            dag->deps[0] = input_dag;
            dag->num_deps = 1;
        }
        dag->G = *D;

        input_dag = dag;

        D++;
    }

    return dag;
}
