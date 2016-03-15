#include "matrix.h"
#include <string.h>
#include <math.h>

void transpose(const int* src, int* dst, int rows_src, int cols_src)
{
  for(int i = 0; i < rows_src; i++)
    for(int j = 0; j < cols_src; j++)
      dst[j*rows_src + i] = src[i*cols_src + j];
}

/*
  Transposes a (m_src x n_src) block of SRC(0,0) into DST(0,0).
  Where src has cols_src columns and dst has cols_dst columns.
*/
void transpose_block(const int* src, int* dst, int m_src, int n_src, int cols_src, int cols_dst)
{
  for(int i = 0; i < m_src; i++)
    for(int j = 0; j < n_src; j++)
      dst[j*cols_dst + i] = src[i*cols_src + j];
}

/*
  Copies a block from src (scols columns)
  starting at (0, 0) and of size (m x n)
  to dst (dcols columns) at (0, 0).
*/
void submat_cpy(const int* src, int scols, int m, int n, int* dst, int dcols)
{
  for(int r = 0; r < m; r++) {
      memcpy(dst, src, sizeof(int)*n);
      src += scols;
      dst += dcols;
  }
}

void rref(int *mat, int rows, int cols, int w)
{
  int i, j, k, x, rs2;
  int row_start, tmp, inverse;

  int r, c;
  r = 0;
  c = 0;

  /* First -- convert into upper triangular  */
  while (r < rows && c < cols)
  {
    row_start = cols*r;

    /*
      Swap rows if we have a zero r,c element.
      If we can't swap, move onto next col.
    */

    if (mat[row_start+c] == 0) {
      for (i = r+1; i < rows && mat[cols*i+c] == 0; i++) ;
      if (i == rows) {
        c++;
        continue;
      }
      rs2 = i*cols;
      for (k = 0; k < cols; k++) {
        tmp = mat[row_start+k];
        mat[row_start+k] = mat[rs2+k];
        mat[rs2+k] = tmp;
      }
    }

    /* Multiply the row by 1/element r,c  */
    tmp = mat[row_start+c];
    if (tmp != 1) {
      inverse = galois_single_divide(1, tmp, w);
      for (j = 0; j < cols; j++) {  //TODO(all): can we start i at c?
        mat[row_start+j] = galois_single_multiply(mat[row_start+j], inverse, w);
      }
    }

    /* Now for each row i > r, subtract A_ic*Ar from Ai  */
    k = row_start+c;
    for (i = r+1; i < rows; i++) {
      k += cols; //mat[k] == mat[i, c]
      if (mat[k] != 0) {
        if (mat[k] == 1) {
          rs2 = cols*i;
          for (x = 0; x < cols; x++) { //TODO(all): x can be started at c?
            mat[rs2+x] ^= mat[row_start+x];
          }
        } else {
          tmp = mat[k];
          rs2 = cols*i;
          for (x = 0; x < cols; x++) { //TODO(all): x can be started at c?
            mat[rs2+x] ^= galois_single_multiply(tmp, mat[row_start+x], w);
          }
        }
      }
    }

    r++;
    c++;
  }

  /*
    Now the matrix is upper triangular.
    Back-substitute.
  */

  for (i = rows-1; i >= 0; i--) {
    row_start = i*cols;

    // find pivot in row i
    int k;
    for (k = i; k < cols && mat[row_start + k] == 0; k++) ;
    if (k == cols)
      continue;

    // substitute into row j
    for (j = 0; j < i; j++) {
      rs2 = j*cols;

      if (mat[rs2+k] != 0) {
        tmp = mat[rs2+k];
        for (x = 0; x < cols; x++) {
          mat[rs2+x] ^= galois_single_multiply(tmp, mat[row_start+x], w);
        }
      }
    }
  }

}

/*
  A version of RREF intended for reducing augmented matrices,
  where only the first rref_cols of the matrix is required to be in RREF form.
  TODO: rref should call this, but keeping separate for now for obvious correctness' sake.
*/
void rref_aug(int *mat, int rows, int cols, int rref_cols, int w)
{
  int i, j, k, x, rs2;
  int row_start, tmp, inverse;

  int r, c;
  r = 0;
  c = 0;

  /* First -- convert into upper triangular  */
  while (r < rows && c < rref_cols)
  {
    row_start = cols*r;

    /*
      Swap rows if we have a zero r,c element.
      If we can't swap, move onto next col.
    */

    if (mat[row_start+c] == 0) {
      for (i = r+1; i < rows && mat[cols*i+c] == 0; i++) ;
      if (i == rows) {
        c++;
        continue;
      }
      rs2 = i*cols;
      for (k = 0; k < cols; k++) {
        tmp = mat[row_start+k];
        mat[row_start+k] = mat[rs2+k];
        mat[rs2+k] = tmp;
      }
    }

    /* Multiply the row by 1/element r,c  */
    tmp = mat[row_start+c];
    if (tmp != 1) {
      inverse = galois_single_divide(1, tmp, w);
      for (j = 0; j < cols; j++) {  //TODO(all): can we start i at c?
        mat[row_start+j] = galois_single_multiply(mat[row_start+j], inverse, w);
      }
    }

    /* Now for each row i > r, subtract A_ic*Ar from Ai  */
    k = row_start+c;
    for (i = r+1; i < rows; i++) {
      k += cols; //mat[k] == mat[i, c]
      if (mat[k] != 0) {
        if (mat[k] == 1) {
          rs2 = cols*i;
          for (x = 0; x < cols; x++) { //TODO(all): x can be started at c?
            mat[rs2+x] ^= mat[row_start+x];
          }
        } else {
          tmp = mat[k];
          rs2 = cols*i;
          for (x = 0; x < cols; x++) { //TODO(all): x can be started at c?
            mat[rs2+x] ^= galois_single_multiply(tmp, mat[row_start+x], w);
          }
        }
      }
    }

    r++;
    c++;
  }

  /*
    Now the matrix is upper triangular.
    Back-substitute.
  */

  for (i = rows-1; i >= 0; i--) {
    row_start = i*cols;

    // find pivot in row i
    int k;
    for (k = i; k < cols && mat[row_start + k] == 0; k++) ;
    if (k == cols || k >= rref_cols)
      continue;

    // substitute into row j
    for (j = 0; j < i; j++) {
      rs2 = j*cols;

      if (mat[rs2+k] != 0) {
        tmp = mat[rs2+k];
        for (x = 0; x < cols; x++) {
          mat[rs2+x] ^= galois_single_multiply(tmp, mat[row_start+x], w);
        }
      }
    }
  }

}

/*
  Column-reduced echelon form.
*/
void cref(int *mat, int rows, int cols, int w)
{
  int* mT = (int*)malloc(sizeof(int) * rows * cols);
  transpose(mat, mT, rows, cols);
  rref(mT, cols, rows, w);
  transpose(mT, mat, cols, rows);
  free(mT);
}

/*
  A version of CREF intended for reducing augmented matrices,
  where only the first rref_rows of the matrix are required to be in CREF form.
*/
void cref_aug(int *mat, int rows, int cols, int cref_rows, int w) {
  int* mT = (int*)malloc(sizeof(int) * rows * cols);
  transpose(mat, mT, rows, cols);
  rref_aug(mT, cols, rows, cref_rows, w);
  transpose(mT, mat, cols, rows);
  free(mT);
}

// G: (nblocks*block_dim) x B
// all_D: concatenated blocks, each block_dim x block_dim
// product <-- diag(all_D) * G
void matrix_left_mult_diag(int* product, int *all_D, int *D_idx, const int nblocks, const int block_dim, int *G, int B, const int w){
    const int n = nblocks;
    const int a = block_dim;

    for (int i = 0; i < n; i++) {
        matrix_multiply(product, all_D + a*a*D_idx[i], G, a, a, a, B, w);
        product += a*B;
        G += a*B;
    }
}

void matrix_multiply(int* product, int *m1, int *m2, int r1, int c1, int r2, int c2, int w)
{
  int i, j, k, l;
  for (i = 0; i < r1*c2; i++) product[i] = 0;

  for (i = 0; i < r1; i++) {
    for (j = 0; j < c2; j++) {
      for (k = 0; k < r2; k++) {
        product[i*c2+j] ^= galois_single_multiply(m1[i*c1+k], m2[k*c2+j], w);
      }
    }
  }
}

// G: rows x B
void matrix_data_product_buff(int* G, int B, int rows,
        char** data_ptrs, char** dst_ptrs_inplace, int data_size,
        int w, int buffsize) {
    char** data_ptrs_curr = (char**)malloc(sizeof(char*) * B);
    char** dst_ptrs_curr = (char**)malloc(sizeof(char*) * rows);

    for(int i = 0; i < B; i++)
        data_ptrs_curr[i] = data_ptrs[i];
    for(int i = 0; i < rows; i++)
        dst_ptrs_curr[i] = dst_ptrs_inplace[i];

    int remainder = data_size;
    while (remainder > 0) {
        int bsize = fmin(buffsize, remainder);
        for (int i = 0; i < rows; i++) {
            int* G_i = G + (B * i);
            row_data_product(G_i, data_ptrs_curr, dst_ptrs_curr[i], B, bsize, w);
            dst_ptrs_curr[i] += bsize;
        }
        for (int i = 0; i < B; i++) {
            data_ptrs_curr[i] += bsize;
        }
        remainder -= bsize;
    }
}

/*
  Multiplies matrix_row by multiple stripes of vector data.
  matrix_row is of length row_len, as is data_ptrs.
  Each data_ptr points to data_size bytes: all the stripes of particular vector-component.
  Puts result in dst.
*/
void row_data_product(int* matrix_row, char** data_ptrs, char* dst, int row_len, int data_size, int w)
{
  //int size = num_stripes * (w/8);
  int size = data_size;
  int k = row_len;

  int init;
  char *dptr, *sptr;
  int i;

  if (w != 1 && w != 8 && w != 16 && w != 32) {
    fprintf(stderr, "ERROR: w is not 1, 8, 16 or 32\n");
    exit(1);
  }

  init = 0;

  dptr = dst;

  /* First copy or xor any data that does not need to be multiplied by a factor */

  for (i = 0; i < k; i++) {
    if (matrix_row[i] == 1) {

      sptr = data_ptrs[i];

      if (init == 0) {
        memcpy(dptr, sptr, size);
        //jerasure_total_memcpy_bytes += size;
        init = 1;
      } else {
        galois_region_xor(sptr, dptr, size);
        //jerasure_total_xor_bytes += size;
      }
    }
  }

  /* Now do the data that needs to be multiplied by a factor */

  for (i = 0; i < k; i++) {
    if (matrix_row[i] != 0 && matrix_row[i] != 1) {

      sptr = data_ptrs[i];

      switch (w) {
        case 8:  galois_w08_region_multiply(sptr, matrix_row[i], size, dptr, init); break;
        case 16: galois_w16_region_multiply(sptr, matrix_row[i], size, dptr, init); break;
        case 32: galois_w32_region_multiply(sptr, matrix_row[i], size, dptr, init); break;
      }
      //jerasure_total_gf_bytes += size;
      init = 1;
    }
  }
}
