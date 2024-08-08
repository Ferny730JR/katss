#ifndef KATSS_BOOTSTRAP_H
#define KATSS_BOOTSTRAP_H

#include <stdbool.h>

typedef enum {
	enrichments,
	ikke,
	counting
} KatssAlgorithm;

struct KatssOptions {
	KatssAlgorithm algo;          /** KatssAlgorithm to bootstrap, e.g. enrichments or counting */
	int            kmer;          /** Length of k to bootstrap */
	int            bs_iters;      /** Number of iterations to bootstrap */
	int            ikke_iters;    /** Number of iterations for ikke, if specified in algo */
	int            sample;        /** Percent to sample. Should be a number between 1-100 */
	bool           probabilistic; /** Use probabilistic method or not */
	int            threads;       /** Number of threads to use */
};

struct KatssBootstrapData {
	uint32_t kmer_hash;
	double mean;
	double stdev;
};

struct KatssBootstrap {
	struct KatssBootstrapData *data;
	uint64_t total;
};

typedef struct KatssOptions KatssOptions;
typedef struct KatssBootstrap KatssBootstrap;
typedef struct KatssBootstrapData KatssBootstrapData;


/**
 * @brief Free the resources associated with bootstrap
 * 
 * @param bootstrap 
 */
void katss_free_bootstrap(KatssBootstrap *bootstrap);


/**
 * @brief Bootstrap a specific file.
 * 
 * Bootstrap is available for k-mer counting, enrichments, and ikke. The
 * `ctrl_file` and `opts` arguments are optional. If not options were passed,
 * the default options is to bootstrap k-mer counts 10 times, while sub-sampling
 * 10% of the file on 5-mers. Modify this by providing custom opts.
 * 
 * @param test_file Test file for enrichments/ikke, or the file for counting
 * @param ctrl_file Optional control file used for enrichments/ikke
 * @param opts      Modify options of bootstrapping
 * @return KatssBootstrap* The bootstrap struct with the information
 */
KatssBootstrap *katss_bootstrap(const char *test_file, const char *ctrl_file, KatssOptions *opts);


/**
 * @brief Initialize the default options for a KatssOptions struct
 * 
 */
void katss_init_default_opts(KatssOptions *opts);

#endif // KATSS_BOOTSTRAP_H
