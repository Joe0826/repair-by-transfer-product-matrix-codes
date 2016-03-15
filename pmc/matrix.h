#ifndef MATRIX_H
#define MATRIX_H

//#include "../jerasure-kmg/include/galois.h"
//#include "../jerasure-kmg/include/jerasure.h"
#include "galois.h"
#include "jerasure.h"
#include <assert.h>

#define invert_matrix(a,b,c,d) assert(jerasure_invert_matrix(a,b,c,d) == 0)
#define mprint jerasure_print_matrix

void transpose(const int* src, int* dst, int rows_src, int cols_src);

/*
  Transposes a (m_src x n_src) block of SRC(0,0) into DST(0,0).
  Where src has cols_src columns and dst has cols_dst columns.
*/
void transpose_block(const int* src, int* dst, int m_src, int n_src, int cols_src, int cols_dst);

/*
  Copies a block from src (scols columns)
  starting at (0, 0) and of size (m x n)
  to dst (dcols columns) at (0, 0).
*/
void submat_cpy(const int* src, int scols, int m, int n, int* dst, int dcols);


void rref(int *mat, int rows, int cols, int w);

/*
  A version of RREF intended for reducing augmented matrices,
  where only the first rref_cols of the matrix is required to be in RREF form.
*/
void rref_aug(int *mat, int rows, int cols, int rref_cols, int w);


/*
  Column-reduced echelon form.
*/
void cref(int *mat, int rows, int cols, int w);

/*
  A version of CREF intended for reducing augmented matrices,
  where only the first rref_rows of the matrix are required to be in CREF form.
*/
void cref_aug(int *mat, int rows, int cols, int cref_rows, int w);


void matrix_left_mult_diag(int* product, int *all_D, int *D_idx, const int nblocks, const int block_dim, int *G, int B, const int w);

void matrix_multiply(int* product, int *m1, int *m2, int r1, int c1, int r2, int c2, int w);

void matrix_data_product_buff(int* G, int row_len, int rows,
        char** data_ptrs, char** dst_ptrs_inplace, int data_size,
        int w, int buffsize);

/*
  Multiplies matrix_row by multiple stripes of vector data.
  matrix_row is of length row_len, as is data_ptrs.
  Each data_ptr points to data_size bytes: all the stripes of particular vector-component.
  Puts result in dst.
*/
void row_data_product(int* matrix_row, char** data_ptrs, char* dst, int row_len, int data_size, int w);

#endif
