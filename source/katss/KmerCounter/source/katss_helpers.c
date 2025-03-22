#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#include "memory_utils.h"
#include "katss.h"

void
katss_init_options(KatssOptions *opts)
{
	opts->kmer = 0;
	opts->iters = 1;
	opts->threads = 1;
	opts->normalize = false;
	opts->sort_enrichments = true;

	opts->bootstrap_iters = 0;
	opts->bootstrap_sample = 25000;

	opts->probs_algo = KATSS_PROBS_NONE;
	opts->probs_ntprec = -1;
	opts->seed = -1;

	opts->enable_warnings = true;
	opts->verbose_output = false;
}

int
katss_parse_options(KatssOptions *opts)
{
	/*================= Catch all possible errors in options =================*/
	/* Check that k-mer is between 1-16 */
	if((opts->kmer < 1 || 16 < opts->kmer) && opts->enable_warnings)
		error_message("KatssOptions: kmer=(%d) must be between 1-16", opts->kmer);
	if(opts->kmer < 1 || 16 < opts->kmer)
		return 1;

	/* Check that iters is within range */
	if(opts->iters < 1 && opts->enable_warnings)
		error_message("KatssOptions: iters=(%d) must be greater than 0", opts->iters);
	if(opts->iters < 1)
		return 1;
	
	/* Check that iters is less than 4^k */
	if(opts->iters > 1ULL << (2*opts->kmer))
		error_message("KatssOptions: iters=(%d) must be less than 4^kmer", opts->iters);
	if(opts->iters > 1ULL << (2*opts->kmer))
		return 1;

	/* Check that threads is greater than 0 */
	if(opts->threads < 0 && opts->enable_warnings)
		error_message("KatssOptions: threads=(%d) must be a non-negative number", opts->threads);
	if(opts->threads < 1)
		return 1;
	
	/* Check that bootstrap iteration is greater than or equal to 0 */
	if(opts->bootstrap_iters < 0 && opts->enable_warnings)
		error_message("KatssOptions: bootstrap_iters=(%d) must be non-negative", opts->bootstrap_iters);
	if(opts->bootstrap_iters < 0)
		return 1;
	
	/* Bootstrap subsample percentage must be within bounds */
	if((opts->bootstrap_sample < 1 || opts->bootstrap_sample > 100000) && opts->enable_warnings)
		error_message("KatssOptions: bootstrap_sample=(%d) must be in range of 1-100000", opts->bootstrap_sample);
	if(opts->bootstrap_sample < 1 || opts->bootstrap_sample > 100000)
		return 1;
	
	/*================= Update values =================*/
	if(opts->probs_ntprec == -1)
		opts->probs_ntprec = (int)round(sqrt((double)opts->kmer));
	if(opts->seed < 0)
		opts->seed = time(NULL);
	
	/*========================== Passed all checks  ==========================*/
	return 0;
}

KatssData *
katss_init_kdata(int kmer)
{
	uint64_t total = 1ULL << (2*kmer);
	KatssData *kdata = s_malloc(sizeof *kdata);
	kdata->num_kmers = total;
	kdata->kmers = s_calloc(total, sizeof *kdata->kmers);

	return kdata;
}

void
katss_free_kdata(KatssData *data)
{
	free(data->kmers);
	free(data);
}
