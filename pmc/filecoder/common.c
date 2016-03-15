#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include "MBR_MSR.h"
#include "common.h"

// Note, MSR the G returned is non-systematic,
// but for MBR it is systematic.
void make_code(CodingTechnique tech, int n, int k, int d, int w,
        int *a, int *B, int** G, int** psi_h,
        int** lambdas_h, int** M_h) {

    switch(tech) {
        case MSR:
            if (d < 2*k - 2) {
                fprintf(stderr,  "  > MSR requires d >= 2k - 2\n");
                exit(0);
            }
            *a   = d-k+1;
            *B   = k*(*a);

            *M_h   = make_m_MSR(k, d);
            *psi_h = talloc(int, n * d);
            *lambdas_h = talloc(int, n);

            make_psi_MSR_ext(*psi_h, *lambdas_h, n, k, d, w); // n x d
            break;
        case MBR:
            *a   = d;
            *B   = k*d - k*(k-1)/2;

            *M_h   = make_m_MBR(k, d);
            *psi_h = make_psi_MBR(n, k, d, w); // n x d
            break;
        default:
            fprintf(stderr,  "  > Unkown coding technique!\n");
            exit(-1);
            break;
    }
    *G   = generator_for_PMC(*psi_h, *M_h, n, d, *a, *B);
    //mprint(*G, n* *a, *B, w);
    //*G_sys_h    = talloc(int, n* (*a) * (*B));
    //memcpy(*G_sys_h, *G_orig_h, sizeof(int) * n*(*a)*(*B));
    //make_systematic(*G_sys_h, n, *a, *B, w);
}

void read_chunk_file(char* chunk, int chunk_size,
        char* code_dir, char* fname, int node, int sym) {
    char chunk_fname[256];
    sprintf(chunk_fname, "%s/%s_node%d_sym%d", code_dir, fname, node, sym);
    dprintf("reading chunk: %s\n", chunk_fname);
    FILE* fchunk = fopen(chunk_fname, "rb");
    if (fchunk == NULL) {
        fprintf(stderr,  "  > Unable to open chunk file: %s\n", chunk_fname);
        exit(1);
    }
    fread(chunk, 1, chunk_size, fchunk);
    fclose(fchunk);
}

int* make_all_mus(const int n, const int k, const int d, const int a,
        const int w, CodingTechnique tech) {
    int* all_mus = talloc(int, n*a);
    for (int i = 0; i < n; i++) {
        make_mu_inplace(all_mus + i*a, i, k, d, w, tech);
    }
    return all_mus;
}

void make_mu_inplace(int* mu,
        const int f, const int k, const int d, const int w,
        const CodingTechnique tech) {
    switch(tech) {
        case MSR:
            make_mu_MSR_inplace(mu, f, k, d, w);
            break;
        case MBR:
            make_mu_MBR_inplace(mu, f, k, d, w);
            break;
    }
}

int* make_mu(
        const int f, const int k, const int d, const int a, const int w,
        CodingTechnique tech) {
    int* mu = talloc(int, a);
    make_mu_inplace(mu, f, k, d, w, tech);
    return mu;
}

void write_meta(FILE *fp,
        char* fname, const int fsize, const int msg_size,
        const int n, const int k, const int d, const int a, const int B, const int w,
        CodingTechnique tech, RbtType rbt) {
    fprintf(fp, "%s\n", fname);
    fprintf(fp, "%d\n", fsize);
    fprintf(fp, "%d\n", msg_size);
    fprintf(fp, "%d\n", n);
    fprintf(fp, "%d\n", k);
    fprintf(fp, "%d\n", d);
    fprintf(fp, "%d\n", a);
    fprintf(fp, "%d\n", B);
    fprintf(fp, "%d\n", w);
    fprintf(fp, "%d\n", tech);
    fprintf(fp, "%d\n", rbt);
}

void read_meta(FILE *fp,
        char* fname, int* fsize, int* msg_size,
        int* n, int* k, int* d, int* a, int* B, int* w,
        CodingTechnique* tech, RbtType* rbt) {
    fscanf(fp, "%s\n", fname);
    fscanf(fp, "%d\n", fsize);     // original file size
    fscanf(fp, "%d\n", msg_size);  // encoded file size (padded)
    fscanf(fp, "%d\n", n);
    fscanf(fp, "%d\n", k);
    fscanf(fp, "%d\n", d);
    fscanf(fp, "%d\n", a);
    fscanf(fp, "%d\n", B);
    fscanf(fp, "%d\n", w);
    fscanf(fp, "%d\n", tech);
    fscanf(fp, "%d\n", rbt);
}

void print_meta( char* fname, int fsize, int msg_size,
        int n, int k, int d, int a, int B, int w,
        CodingTechnique tech, RbtType rbt) {
    printf("n: %d, k: %d, d: %d, w: %d, alpha: %d, B: %d, MBR: %d, RBT: %d\n", n, k, d, w, a, B, tech, rbt);
    printf("file size: %d, msg_size: %d\n", fsize, msg_size);
}

void write_G(FILE *fp, const int* G,
        const int n, const int k, const int d, const int a, const int B, const int w,
        CodingTechnique tech, RbtType rbt) {
    fwrite(&n, sizeof(int), 1, fp);
    fwrite(&k, sizeof(int), 1, fp);
    fwrite(&d, sizeof(int), 1, fp);
    fwrite(&a, sizeof(int), 1, fp);
    fwrite(&B, sizeof(int), 1, fp);
    fwrite(&w, sizeof(int), 1, fp);
    fwrite(&tech, sizeof(CodingTechnique), 1, fp);
    fwrite(&rbt, sizeof(RbtType), 1, fp);

    fwrite(G, sizeof(int), n*a*B, fp);
}

// Returns 1 if fp matches code (read successful),
// else returns 0.
int read_G(FILE *fp, int* G,
        const int n, const int k, const int d, const int a, const int B, const int w,
        CodingTechnique tech, RbtType rbt) {
    int _n,_k,_d,_a,_B,_w;
    CodingTechnique _tech;
    RbtType _rbt;
    fread(&_n, sizeof(int), 1, fp);
    fread(&_k, sizeof(int), 1, fp);
    fread(&_d, sizeof(int), 1, fp);
    fread(&_a, sizeof(int), 1, fp);
    fread(&_B, sizeof(int), 1, fp);
    fread(&_w, sizeof(int), 1, fp);
    fread(&_tech, sizeof(CodingTechnique), 1, fp);
    fread(&_rbt, sizeof(RbtType), 1, fp);
    if (_n == n && _k == k  && _d == d && _a == a && _B == B && _w == w &&
            _tech == tech && _rbt == rbt) {
        fread(G, sizeof(int), n*a*B, fp);
        return 1;
    } else {
        return 0;
    }

}

