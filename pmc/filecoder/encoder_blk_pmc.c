#include <unistd.h>
#include <sys/time.h>
#include <time.h>
#include <sys/stat.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>

#include "MBR_MSR.h"
#include "product_matrix.h"
#include "matrix.h"
#include "dag_helper.h"
#include "rbt.h"
#include "common.h"

int main (int argc, char **argv) {

    printf("WARNING: This is the BLOCKED encoder, which is not yet compatable with collector/repair.\n");
    printf("FOR TESTING USE ONLY.\n");

    FILE *fp;               // file pointers
    CodingTechnique tech;   // coding technique (parameter)
    RbtType rbt = RBT_NONE; // repair-by-transfer type (parameter)
    int n, k,d, w, a, B ;   // PMC parameters
    int use_sparse = 0;

    /** Creation of file name variables. */
    char *curdir, *fname;
    char* code_dir = "coded_files";

    /** Buffers for data and coding matrices. */
    int  *G;
    int *lambdas;
    int *psi;
    int *M;
    int *G_sys = NULL;
    char chunk_fname[256];

    const int total_buffsize = 1 << 18; //16384; //16K
    int buffsize; // will be set after B
    int ioblock; // to be set below.

    /* Repair-by-transfer data. */
    int* all_mus;
    int* all_P;
    int* all_P_inv;
    int* nodes_P_idx;

    /** Timing. */
    double t_IO = 0.0;
    double tv_coding = 0.0; // TODO: transition all timings away from clock()
    double t_overhead = 0.0;
    double t_precode_rbt = 0.0;
    double t_precode_dc = 0.0;
    double t_precode = 0.0;
    double t_postcode_rbt = 0.0;
    double tv_encode_only = 0.0;
    double tv_del = 0.0;
    double t_total = 0.0;
    clock_t t, t2;

	struct timeval tv1, tv2, tv3;

    clock_t start_t = clock();

    /* Error check Arguments*/
    if (argc < 7) {// inputfile n k d
        fprintf(stderr,  "usage: n k d w {MSR/MBR} [-rbt {CYC/SYS}] [-sp] inputfile\n");
        fprintf(stderr,  "\t-sp\t: use sparse precoder\n");
        exit(1);
    }

    /* Conversion of parameters and error checking */
    if (sscanf(argv[1], "%d", &n) == 0 || n <= 0) {
        fprintf(stderr,  "  > Invalid value for n\n");
        exit(1);
    }
    if (sscanf(argv[2], "%d", &k) == 0 || k < 0) {
        fprintf(stderr,  "  > Invalid value for k\n");
        exit(1);
    }
    if (sscanf(argv[3],"%d", &d) == 0 || d <= 0) {
        fprintf(stderr,  "  > Invalid value for d.\n");
        exit(1);
    }
    if (sscanf(argv[4],"%d", &w) == 0 || w <= 0) {
        fprintf(stderr,  "  > Invalid value for w.\n");
        exit(1);
    }

    /* Setting of coding technique and error checking */
    if (strcmp(argv[5], "MSR") == 0) {
        tech = MSR;
    } else if (strcmp(argv[5], "MBR") == 0) {
        tech = MBR;
    } else {
        fprintf(stderr,  "  > Not a valid coding technique. Choose one of the following: {MSR, MBR}.\n");
        exit(1);
    }

    /* Optional repair-by-transfer. */
    int p;
    if (strcmp(argv[6], "-rbt") == 0) {
        if (strcmp(argv[7], "CYC") == 0) {
            rbt = RBT_CYC;
        } else if (strcmp(argv[7], "SYS") == 0) {
            rbt = RBT_SYS;
        } else if (sscanf(argv[7],"SYS%d", &p) == 1 && p >= 0) {
            rbt = RBT_SYSP + p;
            printf("using RBT SYS%d\n", p);
        } else {
            fprintf(stderr,  "  > Not a valid repair-by-transfer type. Choose one of the following: {CYC, SYS}.\n");
            exit(1);
        }
        argv += 2;
    }

    /* Optional flag: use sparse precoder. */
    if (strcmp(argv[6], "-sp") == 0) {
        printf("Using sparse precoder.\n");
        use_sparse = 1;
        argv++;
    } else {
        printf("Using dense precoder.\n");
    }

    /* Make code. */
    t = clock();
    make_code(tech, n, k, d, w,
        &a, &B, &G, &psi, &lambdas, &M);

    if (rbt != RBT_NONE) {
        /* Prepare RBT data. */
        int* all_mus = make_all_mus(n, k, d, a, w, tech);
        gen_rbt_data(&all_P, &all_P_inv, &nodes_P_idx, all_mus, n, a, w, rbt);
    }

    /* Get current working directory for construction of file names */
    curdir = (char*)malloc(sizeof(char)*1000);
    getcwd(curdir, 1000);

    /* Create coded directory */
    t = clock();
    int i = mkdir(code_dir, S_IRWXU);
    if (i == -1 && errno != EEXIST) {
        fprintf(stderr, "  > Unable to create %s directory.\n", code_dir);
        exit(0);
    }
    t_IO += clock()-t;

    /* Open file and error check */
    fp = fopen(argv[6], "rb");
    if (fp == NULL) {
        fprintf(stderr,  "  > Unable to open file.\n");
        exit(0);
    }
    fname = argv[6];

    /* Determine size of file : needs to be a multiple for B*sizeof(long). */
    int fsize, proper_fsize;
    struct stat status;
    stat(argv[6], &status);
    fsize        =  status.st_size;
    proper_fsize = fsize;
    if (fsize % (B*sizeof(long)) != 0)
        proper_fsize += (B*sizeof(long)) - fsize % (B*sizeof(long));
    t_overhead += clock() - t;

    int msg_size = proper_fsize;
    int chunk_size = msg_size / B;

    /* Determine buffsize. */
    buffsize = sizeof(long) * ceil(total_buffsize / (float) (sizeof(long) * B));
    ioblock =  buffsize;
    printf("buffsize: %d\n", buffsize);

    /* Print metadata. */
    print_meta(fname, fsize, msg_size,
        n, k, d, a, B, w, tech, rbt);

    if (use_sparse) {

        /* Read in input file. */
        char** src_ptrs = talloc(char*, B);
        int remainder = fsize;
        t = clock();
        for (int i = 0; i < B; i++){
            int rsize = fmin(chunk_size, remainder);
            aligned_malloc_pmc(src_ptrs + i, chunk_size);
            memset(src_ptrs[i], 0, chunk_size);
            fread(src_ptrs[i], 1, rsize, fp);
            remainder -= rsize;
        }
        fclose(fp);
        t_IO += clock()-t;

        /* Precoding. */

        // Decoding "from" first k nodes.
        t = clock();
        int* sys_nodes = range(0, k);
        GeneratorDAG* dc_precode = NULL;  // Only for MSR.
        if (tech == MSR) {
            // For sparse MSR: First, precode by "data-collecting" from
            // the first-k nodes.
            // Then encode using the (non-systematic, sparse) gen. mtx.
            dc_precode = make_MSR_decoder_dag(psi, lambdas, M, sys_nodes, k, d, w);
        }
        // At this point, G is either non-systematic (if MSR, and sparse precoding)
        // or sparse and systematic (if MBR).


        // Repair-by-transfer precoding:
        // First rbt_decode, then normal data-collect.
        GeneratorDAG* full_precode;
        GeneratorDAG* rbt_precode = NULL;
        if (rbt == RBT_NONE) {
            full_precode = dc_precode;
        } else {
            if (tech == MBR) {
                int* Q = make_rbt_precoder_matrix(all_P, nodes_P_idx, G, k, a, B, w);
                rbt_precode = make_transform_dag(NULL, Q, B, "MBR_precoder");
                full_precode = rbt_precode;
            } else {
                rbt_precode = rbt_decode_dc(all_P_inv, sys_nodes, k, a, w);
                full_precode = concat_dag(rbt_precode, dc_precode);
            }
        }
        t_overhead += clock()-t;

        if (full_precode != NULL) {
            // Perform the precoding.
            char** precoded_ptrs = talloc(char*, B);

            t = clock();
            gettimeofday(&tv1, NULL);
            apply_dag_striped_buff(full_precode, src_ptrs, precoded_ptrs, B, chunk_size, w, buffsize);
            src_ptrs = precoded_ptrs;

            gettimeofday(&tv2, NULL);
            tv_coding += tv2.tv_sec - tv1.tv_sec + (tv2.tv_usec - tv1.tv_usec) / 1000000.0;

            t_precode = clock()-t;

            if (rbt != RBT_NONE) {
                t_precode_rbt = rbt_precode->end_t - rbt_precode->start_t;
                t_precode_dc = full_precode->end_t - rbt_precode->end_t;
            } else {
                t_precode_dc = t_precode;
            }
        }
        //print_dag(full_precode);

        /* Encoding. */
        int* G_i = G;
        char** dst_ptrs = talloc(char*, a);
        char** tmp_ptrs = talloc(char*, a);

        for (int j = 0; j < a; j++) {
            aligned_malloc_pmc(dst_ptrs + j, chunk_size);
            aligned_malloc_pmc(tmp_ptrs + j, chunk_size);
        }

        for (int i = 0; i < n; i++) {
            /* Encode alpha chunks for node i, into dst_ptrs. */

            // First normal encoding (regular generator mtx multiply)
            gettimeofday(&tv1, NULL);
            matrix_data_product_buff(G_i, B, a,
                    src_ptrs, dst_ptrs, chunk_size,
                    w, buffsize);
            gettimeofday(&tv2, NULL);
            tv_del =  tv2.tv_sec - tv1.tv_sec + (tv2.tv_usec - tv1.tv_usec) / 1000000.0;
            tv_coding += tv_del;
            tv_encode_only += tv_del;

            // Then RBT-postcoding (optionally)
            if (rbt != RBT_NONE) {
                char** rbt_src = talloc(char*, a);
                for (int j = 0; j < a; j++) {
                    rbt_src[j] = dst_ptrs[j];
                    dst_ptrs[j] = tmp_ptrs[j];
                }

                int* P = all_P + (a*a)*nodes_P_idx[i];
                t2 = clock();
                gettimeofday(&tv1, NULL);
                matrix_data_product_buff(P, a, a,
                        rbt_src, dst_ptrs, chunk_size,
                        w, buffsize);
                gettimeofday(&tv2, NULL);
                tv_del =  tv2.tv_sec - tv1.tv_sec + (tv2.tv_usec - tv1.tv_usec) / 1000000.0;
                tv_coding += tv_del;
                t_postcode_rbt += clock() - t2;

                for (int j = 0; j < a; j++)
                    tmp_ptrs[j] = rbt_src[j];
            }

            t = clock();
            // Write to alpha files.
            for (int j = 0; j < a; j++) {
                sprintf(chunk_fname, "%s/%s_node%d_sym%d", code_dir, fname, i, j);
                FILE* fchunk = fopen(chunk_fname, "wb");
                fwrite(dst_ptrs[j], 1, chunk_size, fchunk);
                fclose(fchunk);
            }
            t_IO += clock()-t;

            G_i += a*B;
        }
    } else {

        /* Dense encoding. */

        // Try to load a matching precomputed code.
        // (currently only used by the dense strategy).
        int preloadedG = 0;
        FILE* fG = fopen("G.code", "rb");
        if (fG != NULL) {
            preloadedG = read_G(fG, G, n, k, d, a, B, w, tech, rbt);
            fclose(fG);
        }


        if (! preloadedG) {
            // Need to compute G.

            t = clock();
            if (tech == MSR) {
                make_systematic(G, n, a, B, w);
            }

            if (rbt != RBT_NONE) {
                int* rbt_precode;
                int* rbt_postcode;
                if (tech == MBR) {
                    rbt_precode = make_rbt_precoder_matrix(all_P, nodes_P_idx, G, k, a, B, w);
                } else {
                    int* sys_nodes = range(0, k);
                    rbt_precode = rbt_shuffle_mtx(all_P_inv, sys_nodes, k, a, w); // TODO: matrix_right_mult
                }
                int* G_tmp = talloc(int, n*a*B);
                matrix_left_mult_diag(G_tmp, all_P, nodes_P_idx, n, a, G, B, w);
                matrix_multiply(G, G_tmp, rbt_precode, n*a, B, B, B, w);
            }
            t_overhead += clock()-t;

            FILE* fG = fopen("G.code", "wb");
            write_G(fG, G,
                n, k, d, a, B, w, tech, rbt);
            fclose(fG);
        } else {
            printf("Using precomputed G.code\n");
        }
        // G is now the "one-true-G". (systematic, including RBT precode/postcoding).

        // print sparsity of G
        int nsparse = 0;
        int nsparse_sys = 0;
        for (int i=0; i < n*a*B; i++) {
            if (G[i] == 0 || G[i] == 1) {
                nsparse++;
                if (i >= k*a*B)
                    nsparse_sys++;
            }
        }
        float sp_all = (float) nsparse / (n*a*B);
        float sp_sys = (float) nsparse_sys / ((n-k)*a*B);
        //printf("sparsity (all): %f\n", sp_all);
        printf("sparsity (non-sys): %f\n", sp_sys);

        //mprint(G, n*a, B, w);

        /* Read in input file BLOCK. */
        int file_remainder = fsize;
        int chunk_remainder = chunk_size;
        char** src_ptrs = talloc(char*, B);
        for (int i = 0; i < B; i++){
            aligned_malloc_pmc(src_ptrs + i, ioblock);
        }

        char** dst_ptrs = talloc(char*, a);
        for (int i = 0; i < a; i++)
            aligned_malloc_pmc(dst_ptrs + i, ioblock);

        int first_block = 1;
        while (chunk_remainder > 0) {
            t = clock();
            for (int i = 0; i < B; i++) {
                int rsize = fmin(ioblock, file_remainder);
                fread(src_ptrs[i], 1, rsize, fp);
                if (rsize < ioblock)
                    memset(src_ptrs[i] + rsize, 0, ioblock - rsize);
                file_remainder -= rsize;
            }
            chunk_remainder -= ioblock;
            t_IO += clock()-t;


            /* Begin coding.*/
            for (int i = 0; i < n; i++) {
                int* G_i = G + i*a*B;

                gettimeofday(&tv1, NULL);
                matrix_data_product_buff(G_i, B, a,
                        src_ptrs, dst_ptrs, ioblock,
                        w, buffsize);
                gettimeofday(&tv2, NULL);
                tv_del = tv2.tv_sec - tv1.tv_sec + (tv2.tv_usec - tv1.tv_usec) / 1000000.0;
                tv_encode_only += tv_del;
                tv_coding += tv_del;

                // Write to alpha files.
                t = clock();
                for (int j = 0; j < a; j++) {
                    sprintf(chunk_fname, "%s/%s_node%d_sym%d", code_dir, fname, i, j);
                    FILE* fchunk;
                    if (first_block)
                        fchunk = fopen(chunk_fname, "wb");
                    else
                        fchunk = fopen(chunk_fname, "ab");
                    fwrite(dst_ptrs[j], 1, ioblock, fchunk);
                    fclose(fchunk);
                }
                t_IO += clock()-t;
            }
            first_block = 0;

        }
        fclose(fp);

    }

    /* Create metadata file */
    char meta_fname[256];
    sprintf(meta_fname, "%s/%s_metadata", code_dir, fname);
    fp = fopen(meta_fname, "wb");
    write_meta(fp,
        fname, fsize, msg_size,
        n, k, d, a, B, w, tech, rbt);
    fclose(fp);

    /* Free allocated memory */
    free(G);
    free(curdir);

    t_total = clock() - t_total;
    /* Print stats. */
    printf("Bytes read       : %d\n",fsize);
    printf("Bytes written    : %d\n",chunk_size*a*n);
    printf("Overhead [s]     : %.5f\n",t_overhead/CLOCKS_PER_SEC);
    printf("Precode only [s] : %.5f\n",t_precode/CLOCKS_PER_SEC);
    printf("Coding Time [s]  : %.10f\n",tv_coding);
    printf("IO Time [s]      : %.5f\n",t_IO/CLOCKS_PER_SEC);
    printf("Total Time [s]   : %.5f\n",t_total/CLOCKS_PER_SEC);

    printf("\n-- Coding Time Breakdown --\n");
    printf("1. Precode RBT\t: %.5f\n", t_precode_rbt/CLOCKS_PER_SEC);
    printf("2. Precode DC\t: %.5f\n", t_precode_dc/CLOCKS_PER_SEC);
    printf("3. Encode\t: %.5f\n", tv_encode_only);
    printf("4. Postcode RBT\t: %.5f\n", t_postcode_rbt/CLOCKS_PER_SEC);
}
