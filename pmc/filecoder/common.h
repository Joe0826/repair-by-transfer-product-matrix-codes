#ifndef COMMON_H
#define COMMON_H

#include "rbt.h"

#define mprint jerasure_print_matrix
#define talloc(type, num) (type *) malloc(sizeof(type)*(num))

static int verbose = 0;
#ifdef DEBUG
#define dprintf printf
#else
#define dprintf if(NULL) printf
#endif

typedef enum {MSR, MBR} CodingTechnique;

/*
 * Helper function to construct relevant matrices for an (n,k,d,w) code.
 * Populates generator matrix (systematic for MBR, non-systematic for MSR).
 * the encoding matrix psi, lambdas (for MSR), and the
 * structure-specifying message matrix M.
 */
void make_code(CodingTechnique tech, int n, int k, int d, int w,
        int *a, int *B, int** G, int** psi_h,
        int** lambdas_h, int** M_h);


void read_chunk_file(char* chunk, int chunk_size,
        char* code_dir, char* fname, int node, int sym);

int* make_all_mus(const int n, const int k, const int d, const int a,
        const int w, CodingTechnique tech);

void make_mu_inplace(int* mu,
        const int f, const int k, const int d, const int w,
        const CodingTechnique tech);

int* make_mu(
        const int f, const int k, const int d, const int a, const int w,
        CodingTechnique tech);

void write_meta(FILE *fp,
        char* fname, const int fsize, const int msg_size,
        const int n, const int k, const int d, const int a, const int B, const int w,
        CodingTechnique tech, RbtType rbt);

void read_meta(FILE *fp,
        char* fname, int* fsize, int* msg_size,
        int* n, int* k, int* d, int* a, int* B, int* w,
        CodingTechnique* tech, RbtType* rbt);

void print_meta( char* fname, int fsize, int msg_size,
        int n, int k, int d, int a, int B, int w,
        CodingTechnique tech, RbtType rbt);

void write_G(FILE *fp, const int* G,
        const int n, const int k, const int d, const int a, const int B, const int w,
        CodingTechnique tech, RbtType rbt);

int read_G(FILE *fp, int* G,
        const int n, const int k, const int d, const int a, const int B, const int w,
        CodingTechnique tech, RbtType rbt);
#endif
