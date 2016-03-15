#include <unistd.h>
#include <sys/time.h>
#include <time.h>
#include <sys/stat.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <unistd.h>

#include "MBR_MSR.h"
#include "product_matrix.h"
#include "rbt.h"
#include "common.h"

int main (int argc, char **argv) {
    FILE *fp;
    CodingTechnique tech;   // coding technique (parameter)
    RbtType rbt;            // repair-by-transfer type (parameter)
    int n, k, d, w, a, B;   // PMC parameters
    int f;                  // the index of the failed node (to help/fix).
    int msg_size;  // original file size
    int fsize;     // encoded file size (padded)

    /** Creation of file name variables. */
    char fname[256];
    char chunk_fname[256];
    char* code_dir = "coded_files";
    char* helper_dir = "helper_files";

    printf("usage: iotest [metadata fname] {link/read}\n");

    /* Open encoder metadata file and parse. */
    fp = fopen(argv[1], "r");
    if (fp == NULL) {
        fprintf(stderr,  "  > Unable to open file.\n");
        exit(1);
    }
    read_meta(fp,
        fname, &fsize, &msg_size,
        &n, &k, &d, &a, &B, &w, &tech, &rbt);
    fclose(fp);

    /* Print metadata. */
    print_meta(fname, fsize, msg_size,
        n, k, d, a, B, w, tech, rbt);

    int chunk_size = msg_size / B;

    /* Create helper directory */
    //int i = mkdir(helper_dir, S_IRWXU);
    //if (i == -1 && errno != EEXIST) {
        //fprintf(stderr, "  > Unable to create %s directory.\n", helper_dir);
        //exit(1);
    //}

    int h_id = 0; // helper-id
    f = 1; //failed-id
    if (strcmp(argv[2], "link") == 0) {
        printf("*** symlinking REPAIR BY TRANSFER chunk from node %d to failed node %d...\n", h_id, f);

        int idx = 0; //get_repair_chunk_idx(h_id, f, n, a, rbt);

        sprintf(chunk_fname, "%s/%s_helper_chunk_from%d_to%d", helper_dir, fname, h_id, f);
        unlink(chunk_fname); // rm existing link, if any
        // to avoid the following DANGEROUS situation:
        // helper_chunk --> coded_file
        // write to helper_chunk... actually writes to coded_file!
        char stored_chunk_fname[256];
        // we need the SOURCE path (in coded_files) to be relative to the TARGET path (helper_files)
        // hence the leading '../'
        sprintf(stored_chunk_fname, "../%s/%s_node%d_sym%d", code_dir, fname, h_id, idx);
        int err = symlink(stored_chunk_fname, chunk_fname);
        if (err != 0) {
            fprintf(stderr,  "  > Unable to symlink helper chunk. ERRNO: %d", errno);
            exit(1);
        }
    } else {
        // Need to compute the repair chunk.
        printf("*** reading all data for repair chunk from node %d to failed node %d...\n", h_id, f);

        char** chunk_ptrs = talloc(char*, a);
        for (int i = 0; i < a; i++)
            chunk_ptrs[i] = talloc(char, chunk_size);
        //[> Allocate for helper-chunk. <]
        char* helper_chunk = talloc(char, chunk_size);

        //[> Read in all stored chunks. <]
        for (int i = 0; i < a; i++) {
            read_chunk_file(chunk_ptrs[i], chunk_size,
                    code_dir, fname, h_id, i);
        }
    }

}
