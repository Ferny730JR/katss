#ifndef KATSS_STRUCTURE_H
#define KATSS_STRUCTURE_H

#include <stdbool.h>
#include "bpp_tables.h"

struct BppOptions {
	int kmer;
	int threads;
	bool seq_windows;
	int window_size;
};

typedef struct BppOptions BppOptions;

void
katss_bpp_init_default_opts(BppOptions *opts);


/**
 * @brief Get the k-mer base-pair probabilities from the test file, normalized
 * by the control file. Set the k-mer length, number of threads, and algorithm
 * in the BppOptions struct.
 * 
 * @param test_file 
 * @param ctrl_file 
 * @param opts 
 * @return kmerHashTable* 
 */
kmerHashTable *
katss_bpp(const char *test_file, const char *ctrl_file, BppOptions *opts);

void
katss_free_bpp(kmerHashTable *bpp);

#endif
