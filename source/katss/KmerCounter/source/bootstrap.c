#include <stdbool.h>
#include <string.h>
#include <math.h>

#include "katss_core.h"
#include "counter.h"
#include "hash_functions.h"
#include "memory_utils.h"
#include "enrichments.h"
#include "bootstrap.h"

void
katss_init_default_opts(KatssOptions *opts)
{
	opts->algo          = counting;
	opts->kmer          = 5;
	opts->bs_iters      = 10;
	opts->ikke_iters    = 1;
	opts->sample        = 10;
	opts->probabilistic = false;
	opts->threads       = 8;
}

static KatssBootstrap *
katss_init_bootstrap(int kmer)
{
	uint64_t total = ((uint64_t)1) << (2 * kmer);
	KatssBootstrap *bootstrap = s_malloc(sizeof *bootstrap);
	bootstrap->data = s_calloc(total, sizeof *bootstrap->data); // init values to 0
	bootstrap->total = total;

	return bootstrap;
}

static void
running_stdev(double value, double *mean, double *stdev, int run)
{
	double tmp_mean = *mean;
	*mean  += (value - tmp_mean) / run;
	*stdev += (value - tmp_mean) * (value - *mean);
}

static double
predict_kmer(char *kseq, KatssCounter *monomer_counts, KatssCounter *dimer_counts)
{
	char monoseq[2] = { 0 };
	char diseq[3]   = { 0 };

	double monoprob = 1;
	double diprob = 1;

	int kmer = strlen(kseq);

	/* Get the cumulative probabilities for overlapping monomers */
	for(int i=1; i<kmer-1; i++) {
		monoseq[0] = kseq[i];
		double count;
		katss_get(monomer_counts, KATSS_DOUBLE, &count, monoseq);
		monoprob *= count/monomer_counts->total;
	}

	/* Get probabilities for all di-mers in k-mer */
	for(int i=0; i<kmer-1; i++) {
		diseq[0] = kseq[i]; diseq[1] = kseq[i+1];
		double count;
		katss_get(dimer_counts, KATSS_DOUBLE, &count, diseq);
		diprob *= count/dimer_counts->total;
	}

	/* Predicted k-mer probability is dinucleotides / overlapping monomers */
	return diprob/monoprob;
}

static int
process_count(const char *file, KatssBootstrap *bootstrap, KatssOptions *opts, int run)
{
	unsigned int kmer = (unsigned int)opts->kmer;
	int smp = opts->sample;
	int thr = opts->threads;

	KatssCounter *counter = katss_count_kmers_bootstrap_mt(file, kmer, smp, thr);
	if(counter == NULL)
		return 1;

	/* Update standard deviation */
	uint64_t count;
	for(uint64_t i=0; i<bootstrap->total; i++) {
		katss_get_from_hash(counter, KATSS_UINT64, &count, (uint32_t)i);
		running_stdev(count, &bootstrap->data[i].mean, &bootstrap->data[i].stdev, run);
		bootstrap->data[i].kmer_hash = (uint32_t)i;
	}

	/* No error encountered */
	katss_free_counter(counter);
	return 0;
}

static int
process_enrichments_prob(const char *test_file, KatssBootstrap *bootstrap, KatssOptions *opts, int run)
{
	int return_code = 1; // assume error in counting
	KatssCounter *test_counts = katss_count_kmers_bootstrap_mt(test_file, opts->kmer, opts->sample, opts->threads);
	if(test_counts == NULL)
		goto exit;
	KatssCounter *mono_counts = katss_count_kmers_bootstrap_mt(test_file, 1, opts->sample, opts->threads);
	if(mono_counts == NULL)
		goto undo_test;
	KatssCounter *dint_counts = katss_count_kmers_bootstrap_mt(test_file, 2, opts->sample, opts->threads);
	if(dint_counts == NULL)
		goto undo_mono;

	char kseq[17];
	return_code = 2; // error in hash table
	for(uint64_t i=0; i<bootstrap->total; i++) {
		double test_frq, ctrl_frq;
		if(katss_get_from_hash(test_counts, KATSS_DOUBLE, &test_frq, (uint32_t)i) != 0)
			goto undo_dint;
		katss_unhash(kseq, (uint32_t)i, opts->kmer, true);
		test_frq = test_frq / test_counts->total;
		ctrl_frq = predict_kmer(kseq, mono_counts, dint_counts);

		if(test_frq == 0.0 || ctrl_frq == 0.0)
			continue;
		
		double rval = test_frq / ctrl_frq;
		running_stdev(rval, &bootstrap->data[i].mean, &bootstrap->data[i].stdev, run);
		bootstrap->data[i].kmer_hash = (uint32_t)i;
	}
	return_code = 0;

undo_dint:
	katss_free_counter(dint_counts);
undo_mono:
	katss_free_counter(mono_counts);
undo_test:
	katss_free_counter(test_counts);
exit:
	return return_code;
}

static int
process_enrichments(const char *test_file, const char *ctrl_file, 
                    KatssBootstrap *bootstrap, KatssOptions *opts, int run)
{
	if(opts->probabilistic)
		return process_enrichments_prob(test_file, bootstrap, opts, run);

	int return_code = 1;
	KatssCounter *test_counts = katss_count_kmers_bootstrap_mt(test_file, opts->kmer, opts->sample, opts->threads);
	if(test_counts == NULL)
		goto exit;
	KatssCounter *ctrl_counts = katss_count_kmers_bootstrap_mt(ctrl_file, opts->kmer, opts->sample, opts->threads);
	if(ctrl_counts == NULL)
		goto undo_test;

	/* Compute enrichments */
	for(uint64_t i=0; i<bootstrap->total; i++) {
		double test_frq, ctrl_frq;
		if(katss_get_from_hash(test_counts, KATSS_DOUBLE, &test_frq, (uint32_t)i) != 0)
			goto undo_ctrl;
		if(katss_get_from_hash(ctrl_counts, KATSS_DOUBLE, &ctrl_frq, (uint32_t)i) != 0)
			goto undo_ctrl;
		if(test_frq == 0.0 || ctrl_frq == 0.0)
			continue;

		test_frq = test_frq / test_counts->total;
		ctrl_frq = ctrl_frq / ctrl_counts->total;

		double rval = test_frq / ctrl_frq;
		running_stdev(rval, &bootstrap->data[i].mean, &bootstrap->data[i].stdev, run);
		bootstrap->data[i].kmer_hash = (uint32_t)i;
	}
	return_code = 0;

undo_ctrl:
	katss_free_counter(ctrl_counts);
undo_test:
	katss_free_counter(test_counts);
exit:
	return return_code;
}

static int
process_ikke_prob(const char *test_file, KatssBootstrap *bootstrap, KatssOptions *opts, int run)
{
	return 1;
}

static int
process_ikke(const char *test_file, const char *ctrl_file, 
             KatssBootstrap *bootstrap, KatssOptions *opts, int run)
{
	if(opts->probabilistic)
		return process_ikke_prob(test_file, bootstrap, opts, run);
	return 1;
}

static int
process_bootstrap_iteration(
	const char *test_file, const char *ctrl_file, KatssBootstrap *bootstrap, KatssOptions *opts, int run)
{
	switch(opts->algo) {
	case counting:    return process_count(test_file, bootstrap, opts, run);
	case enrichments: return process_enrichments(test_file, ctrl_file, bootstrap, opts, run);
	case ikke:        return process_ikke(test_file, ctrl_file, bootstrap, opts, run);
	default: return 1;
	}
}

void
katss_free_bootstrap(KatssBootstrap *bootstrap)
{
	free(bootstrap->data);
	free(bootstrap);
}

static int
bootstrap_compare(const void *a, const void *b)
{
	const KatssBootstrapData *p1 = (KatssBootstrapData *)a;
	const KatssBootstrapData *p2 = (KatssBootstrapData *)b;

	if(p1->mean < p2->mean)
		return 1;
	else if(p1->mean > p2->mean)
		return -1;
	return 0;
}

KatssBootstrap *
katss_bootstrap(const char *test_file, const char *ctrl_file, KatssOptions *opts)
{
	/* Check if options were passed */
	bool options_passed = true;
	if(opts == NULL) {
		options_passed = false;
		opts = s_malloc(sizeof *opts);
		katss_init_default_opts(opts);
	}

	/* Check parameters */
	if(test_file == NULL)
		return NULL;
	if(ctrl_file == NULL && !opts->probabilistic && opts->algo != counting)
		return NULL;
	if(opts->kmer < 1 || 16 < opts->kmer)
		return NULL;
	if(opts->bs_iters < 1 || opts->sample < 1 || 100 < opts->sample)
		return NULL;
	if(opts->algo == ikke && opts->ikke_iters < 1)
		return NULL;

	/* Initialize bootstrap */
	KatssBootstrap *bootstrap = katss_init_bootstrap(opts->kmer);

	/* Process the iterations */
	for(int i=1; i <= opts->bs_iters; i++) {
		if(process_bootstrap_iteration(test_file, ctrl_file, bootstrap, opts, i) != 0) {
			katss_free_bootstrap(bootstrap);
			bootstrap = NULL;
			break;
		}
	}

	/* Finish computing the standard deviation */
	for(uint64_t i = 0; bootstrap != NULL && i < bootstrap->total; i++) {
		bootstrap->data[i].stdev = sqrt(bootstrap->data[i].stdev / (opts->bs_iters - 1));
	}

	if(bootstrap != NULL)
		qsort(bootstrap->data, bootstrap->total, sizeof *bootstrap->data, bootstrap_compare);

	if(!options_passed)
		free(opts);
	return bootstrap;
}
