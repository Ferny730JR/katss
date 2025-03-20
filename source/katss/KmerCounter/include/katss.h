#ifndef KATSS_H
#define KATSS_H

#include <stdbool.h>
#include <stdint.h>

/*==============================================================================
 STRUCT DEFINITIONS
==============================================================================*/

/**
 * @brief Information stored for a specific k-mer
 */
struct KatssDataEntry {
	uint32_t kmer;
	union {
		float rval;
		uint32_t count;
	};
	double pval;
	float stdev;
};
typedef struct KatssDataEntry KatssDataEntry;

/**
 * @brief All k-mer entries
 */
struct KatssData {	
	KatssDataEntry *kmers;
	uint64_t num_kmers;
};
typedef struct KatssData KatssData;

/**
 * @brief Which probabilistic algorithm to use, if any
 */
typedef enum {
	KATSS_PROBS_NONE,
	KATSS_PROBS_REGULAR,
	KATSS_PROBS_USHUFFLE,
	KATSS_PROBS_BOTH
} KatssProbsAlgo;

/**
 * @brief Options used to modify the output of the katss_* functions
 */
struct KatssOptions {
	int kmer;              /** Set length of k-mers */
	int iters;             /** Set the number of iterations for ikke */
	int threads;           /** Set the number of threads to use */
	int normalize;         /** Get the log2 of the enrichments */
	bool sort_enrichments; /** Sort the enrichment based on rval. This will
	                           result data->kmers[0] to be the kmer with the
	                           highest rval, and the following kmers in
	                           decreasing orders based on rval */

	/* bootstrap options */
	int bootstrap_iters;   /** Number of iterations to bootstrap. 0 to not bootstrap */
	int bootstrap_sample;  /** Percent to sample. Should be a number between 
	                           1-100000. Every number represents 0.001%, e.g.,
	                           250000 -> 25.000%, 12345 -> 12.345% */
	
	/* Probabilistic Options */
	KatssProbsAlgo probs_algo;   /* Specify which probabilistic method to use */
	int            probs_ntprec; /* Precision in kmer prediction. Set it as -1 for recommended value */
	int            seed;         /* Seed to use for which random sequences to sample */

	/* Function information */
	bool enable_warnings;        /* Display warnings regarding options */
	bool verbose_output;         /* Display verbose output of calculations */
};
typedef struct KatssOptions KatssOptions;


/*==============================================================================
 FUNCTION DEFINITIONS
==============================================================================*/

/**
 * @brief Initialize all parameters to the default configuration.
 * 
 * @param opts Pointer to KatssOptions struct
 */
void
katss_init_options(KatssOptions *opts);


/**
 * @brief Count all kmers in a dataset
 * 
 * @param path File path to the dataset that will be counter. Can be either a
 * fasta, fastq, or reads file format. Work with files that are gzip compressed.
 * @param opts Options struct to modify the counting algorithm
 * @return KatssData* 
 */
KatssData *
katss_count(const char *path, KatssOptions *opts);


/**
 * @brief Compute the most enriched k-mer in a sequence dataset
 * 
 * @param test File path to the test dataset. Can be either a fasta, fastq, or
 * reads file format. Works with files using gzip compression.
 * @param ctrl File path to the control dataset. This parameter is optional; if
 * there is no control dataset, pass a NULL pointer.
 * @param opts Options struct to modify the enrichment algorithm.
 * @return KatssData* Pointer to computed enrichment values
 */
KatssData *
katss_enrichment(const char *test, const char *ctrl, KatssOptions *opts);


/**
 * @brief Compute the iterative k-mer knockout enrichments
 * 
 * @param test File path to the test dataset
 * @param ctrl File path to the control dataset. This parameter is optional (if
 * using a probabilistic algorithm) - pass a NULL pointer
 * @param opts Options struct to modify the output
 * @return KatssData* Pointer to computed enrichment values
 */
KatssData *
katss_ikke(const char *test, const char *ctrl, KatssOptions *opts);


/**
 * @brief Free all allocated resources to kdata
 * 
 * @param data 
 */
void
katss_free_kdata(KatssData *data);

#endif // KATSS_H
