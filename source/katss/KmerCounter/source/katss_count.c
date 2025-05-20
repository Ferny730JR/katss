#include <stddef.h>
#include <stdlib.h>
#include <math.h>

#include "katss.h"
#include "katss_helpers.h"
#include "memory_utils.h"

#include "counter.h"
#include "ushuffle.h"

static int
compare1(const void *a, const void *b)
{
	const KatssDataEntry *p1 = (KatssDataEntry *)a;
	const KatssDataEntry *p2 = (KatssDataEntry *)b;

	if(p1->count < p2->count) {
		return 1;
	} else if(p1->count > p2->count) {
		return -1;
	} else {
		return 0;
	}
}

static int
compare2(const void *a, const void *b)
{
	const KatssDataEntry *p1 = (KatssDataEntry *)a;
	const KatssDataEntry *p2 = (KatssDataEntry *)b;

	if(isnan(p1->rval) && isnan(p2->rval)) {
		return 0;
	} else if(isnan(p1->rval)) {
		return 1;
	} else if(isnan(p2->rval)) {
		return -1;
	}

	if(p1->rval < p2->rval) {
		return 1;
	} else if(p1->rval > p2->rval) {
		return -1;
	} else {
		return 0;
	}
}

static void
running_stdev(float value, float *mean, float *stdev, int run)
{
	float tmp_mean = *mean;
	*mean  += (value - tmp_mean) / run;
	*stdev += (value - tmp_mean) * (value - *mean);
}

static KatssData *
regular(const char *path, KatssOptions *opts)
{
	/* Compute counts */
	KatssCounter *ctr = katss_count_kmers_mt(path, opts->kmer, opts->threads);
	if(ctr == NULL)
		return NULL;
	
	/* Move counts to KatssData */
	KatssData *counts = katss_init_kdata(opts->kmer);
	for(uint64_t i=0; i<counts->num_kmers; i++) {
		counts->kmers[i].kmer = (uint32_t)i;
		katss_get_from_hash(ctr, KATSS_UINT32, &counts->kmers[i].count, (uint32_t)i);
	}

	/* Free data */
	katss_free_counter(ctr);

	return counts;
}

static KatssData *
ushuffle(const char *path, KatssOptions *opts)
{
	/* Compute shuffled counts */
	KatssCounter *ctr = katss_count_kmers_ushuffle(path, opts->kmer, opts->probs_ntprec);
	if(ctr == NULL)
		return NULL;
	
	/* Move counts to KatssData */
	KatssData *counts = katss_init_kdata(opts->kmer);
	for(uint64_t i=0; i<counts->num_kmers; i++) {
		counts->kmers[i].kmer = (uint32_t)i;
		katss_get_from_hash(ctr, KATSS_UINT32, &counts->kmers[i].count, (uint32_t)i);
	}

	/* Free data */
	free(ctr);

	return counts;
}

static KatssData *
bootstrap_regular(const char *path, KatssOptions *opts)
{
	KatssData *counts = katss_init_kdata(opts->kmer);
	if(counts == NULL)
		return NULL;
	KatssCounter *ctr = NULL;
	unsigned int seed = opts->seed;
	unsigned int kmer = opts->kmer;
	int sample        = opts->bootstrap_sample;
	int threads       = opts->threads;

	/* Process n number of iterations */
	for(int i=1; i<=opts->bootstrap_iters; i++) {
		/* Compute counts */
		ctr = katss_count_kmers_bootstrap_mt(path, kmer, sample, &seed, threads);
		if(ctr == NULL) {
			error_message("katss_count: Failed to get counts on iteration=(%d)", i);
			katss_free_kdata(counts);
			return NULL;
		}

		/* Move counts to KatssData */
		float count;
		for(uint64_t n=0; n<counts->num_kmers; n++) {
			counts->kmers[n].kmer = (uint32_t)n;
			if(katss_get_from_hash(ctr, KATSS_FLOAT, &count, (uint32_t)n) != 0)
				goto error;
			running_stdev(count, &counts->kmers[n].rval, &counts->kmers[n].stdev, i);
		}

		/* Free data */
		katss_free_counter(ctr);
	}

	/* Finish computing stdev */
	if(opts->bootstrap_iters > 1)
		for(uint64_t i=0; i<counts->num_kmers; i++)
			counts->kmers[i].stdev = sqrt(counts->kmers[i].stdev / (opts->bootstrap_iters - 1));

	return counts;

error:
	katss_free_counter(ctr);
	katss_free_kdata(counts);
	return NULL;
}

static KatssData *
bootstrap_ushuffle(const char *path, KatssOptions *opts)
{
	KatssData *counts = katss_init_kdata(opts->kmer);
	if(counts == NULL)
		return NULL;
	KatssCounter *ctr;
	unsigned int seed = opts->seed;
	unsigned int kmer = opts->kmer;
	int klet          = opts->probs_ntprec;
	int sample        = opts->bootstrap_sample;

	/* Process n number of iterations */
	for(int i=1; i<=opts->bootstrap_iters; i++) {
		/* Compute counts */
		ctr = katss_count_kmers_ushuffle_bootstrap(path, kmer, klet, sample, &seed);
		if(ctr == NULL) {
			error_message("katss_count: Failed to get counts on iteration=(%d)", i);
			katss_free_kdata(counts);
			return NULL;
		}
		
		/* Move counts to KatssData */
		float count;
		for(uint64_t n=0; n<counts->num_kmers; n++) {
			counts->kmers[n].kmer = (uint32_t)n;
			katss_get_from_hash(ctr, KATSS_FLOAT, &count, (uint32_t)n);
			running_stdev(count, &counts->kmers[n].rval, &counts->kmers[n].stdev, i);
		}

		/* Free data */
		katss_free_counter(ctr);
	}

	/* Finish computing stdev */
	for(uint64_t n=0; n<counts->num_kmers; n++) {
		counts->kmers[n].stdev = sqrt(counts->kmers[n].stdev / (opts->bootstrap_iters - 1));
	}

	return counts;
}

KatssData *
katss_count(const char *path, KatssOptions *opts)
{
	/* Make sure test file was passed */
	if(path == NULL)
		return NULL;

	/* Parse the options */
	if(katss_parse_options(opts) != 0)
		return NULL;


	/* BEGIN COMPUTATION: No bootstrap */
	KatssData *data = NULL;
	if(opts->bootstrap_iters == 0) {
		switch(opts->probs_algo) {
		case KATSS_PROBS_NONE:
			data = regular(path, opts);  // Just get the counts
			break;

		case KATSS_PROBS_REGULAR:
			if(opts->enable_warnings)
				error_message("katss_count: KATSS_PROBS_REGULAR is not supported");
			return NULL; // can't compute

		case KATSS_PROBS_USHUFFLE:
			data = ushuffle(path, opts); // Get counts of shuffled seq
			break;

		case KATSS_PROBS_BOTH:
			if(opts->enable_warnings)
				error_message("katss_count: KATSS_PROBS_BOTH is not supported");
			return NULL; // can't compute

		default: return NULL; // silence compiler warnings 
		}

	/* BEGIN COMPUTATION: bootstrap */
	} else {
		switch(opts->probs_algo) {
		case KATSS_PROBS_NONE:
			data = bootstrap_regular(path, opts); // bootstrap counts
			break;

		case KATSS_PROBS_REGULAR:
			if(opts->enable_warnings)
				error_message("katss_count: KATSS_PROBS_REGULAR is not supported");
			return NULL; // can't compute

		case KATSS_PROBS_USHUFFLE:
			data = bootstrap_ushuffle(path, opts); // bootstrap ushuffle counts
			break;

		case KATSS_PROBS_BOTH:
			if(opts->enable_warnings)
				error_message("katss_count: KATSS_PROBS_BOTH is not supported");
			return NULL;

		default: return NULL;
		}
	}

	/* If data is NULL, reeturn NULL */
	if(data == NULL)
		return NULL;

	/* Sort if necessary */
	if(opts->sort_enrichments && opts->bootstrap_iters == 0) {
		qsort(data->kmers, data->num_kmers, sizeof *data->kmers, compare1);
	} else if(opts->sort_enrichments) {
		qsort(data->kmers, data->num_kmers, sizeof *data->kmers, compare2);
	}

	/* DONE: return data */
	return data;
}
