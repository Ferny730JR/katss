#include <stddef.h>
#include <stdlib.h>
#include <math.h>

#include "katss.h"
#include "katss_core.h"
#include "katss_helpers.h"
#include "memory_utils.h"

#include "enrichments.h"
#include "t_test.h"

static int
compare(const void *a, const void *b)
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
running_stdev(double value, double *mean, double *stdev, int run)
{
	double tmp_mean = *mean;
	*mean  += (value - tmp_mean) / run;
	*stdev += (value - tmp_mean) * (value - *mean);
}

/**
 * @brief Compute the enrichments of all kmers
 * 
 * @param test Test file to be used for computation
 * @param ctrl Control file to be used for computation
 * @param opts Options to modify output
 * @return KatssData* Data containing rval's
 */
static KatssData *
regular(const char *test, const char *ctrl, KatssOptions *opts)
{
	KatssEnrichments *enr;
	KatssData *enrichments;

	/* Compute enrichments */
	enr = katss_enrichments(test, ctrl, opts->kmer, opts->normalize);
	if(enr == NULL)
		return NULL;

	/* Move enrichments to KatssData */
	enrichments = katss_init_kdata(opts->kmer);
	for(uint64_t i=0; i<enrichments->num_kmers; i++) {
		enrichments->kmers[i].kmer = enr->enrichments[i].key;
		enrichments->kmers[i].rval = (float)enr->enrichments[i].enrichment;
	}

	/* Free data */
	katss_free_enrichments(enr);

	return enrichments;
}

/**
 * @brief Compute the enrichments using the probabilistic method.
 * 
 * @param test Testfile to be used for computation
 * @param opts Options to modify output
 * @return KatssData* Data containing rval's
 */
static KatssData *
probs(const char *test, KatssOptions *opts)
{
	KatssEnrichments *enr = NULL;

	/* Compute probabilistic enrichments */
	enr = katss_prob_enrichments(test, opts->kmer, opts->normalize);
	if(enr == NULL)
		return NULL;

	/* Initialize kdata */
	KatssData *data = katss_init_kdata(opts->kmer);
	if(data == NULL)
		goto exit;
	
	/* Move enrichments to KatssData */
	for(uint64_t i=0; i<data->num_kmers; i++) {
		data->kmers[i].kmer = enr->enrichments[i].key;
		data->kmers[i].rval = (float)enr->enrichments[i].enrichment;
	}

exit:
	katss_free_enrichments(enr);
	return data;
}

/**
 * @brief Compute the enrichments of the shuffled sequences.
 * 
 * @param test Testfile to be used for computation
 * @param opts Options to modify output
 * @return KatssData* Data containing rval's
 */
static KatssData *
ushuffle(const char *test, KatssOptions *opts)
{
	KatssData *data = NULL;
	KatssEnrichments *enr = NULL;
	unsigned int kmer = opts->kmer;
	int klet          = opts->probs_ntprec;
	bool normalize    = opts->normalize;

	/* Compute the counts */
	KatssCounter *test_counts = katss_count_kmers(test, kmer);
	if(test_counts == NULL)
		goto exit_error;
	KatssCounter *shuf_counts = katss_count_kmers_ushuffle(test, kmer, klet);
	if(shuf_counts == NULL)
		goto exit_error;

	/* Compute the enrichments */
	enr = katss_compute_enrichments(test_counts, shuf_counts, normalize);
	if(enr == NULL)
		goto exit_error;
	
	/* Free counts */
	katss_free_counter(test_counts);
	katss_free_counter(shuf_counts);

	/* Initialize kdata */
	data = katss_init_kdata(kmer);
	if(data == NULL)
		goto undo_enrichment;

	/* Move enrichments to KatssData */
	for(uint64_t i=0; i<data->num_kmers; i++) {
		data->kmers[i].kmer = enr->enrichments[i].key;
		data->kmers[i].rval = (float)enr->enrichments[i].enrichment;
	}

	/* Success: return */
	katss_free_enrichments(enr);
	return data;

/* ERRORS ENCOUNTERED */
undo_enrichment:
	katss_free_enrichments(enr);
exit_error:
	katss_free_counter(test_counts);
	katss_free_counter(shuf_counts);
	return NULL;
}

/**
 * @brief Compute the enrichments using both the shuffled and probabilistic
 * method
 * 
 * @param test Testfile to be used for computation
 * @param opts Options to modify output
 * @return KatssData* Data containing rval's
 */
static KatssData *
both(const char *test, KatssOptions *opts)
{
	KatssData *data = NULL;
	KatssEnrichments *shuf = NULL;
	KatssEnrichments *prob = NULL;
	unsigned int kmer = opts->kmer;
	int klet          = opts->probs_ntprec;

	/* Compute the counts */	
	KatssCounter *test_counts = katss_count_kmers_ushuffle(test, kmer, klet);
	KatssCounter *mono_counts = katss_count_kmers_ushuffle(test, 1, klet);
	KatssCounter *dint_counts = katss_count_kmers_ushuffle(test, 2, klet);
	if(test_counts == NULL || mono_counts == NULL || dint_counts == NULL)
		goto exit_error;

	/* Compute the enrichments */
	shuf = katss_compute_prob_enrichments(test_counts, mono_counts, dint_counts, false);
	if(shuf == NULL)
		goto exit_error;

	/* Free counters */
	katss_free_counter(test_counts);
	katss_free_counter(mono_counts);
	katss_free_counter(dint_counts);

	/* Compute the probabilistic enrichments */
	prob = katss_prob_enrichments(test, kmer, false);
	if(prob == NULL)
		goto exit;

	/* Compute rval from both probabilistic methods */
	data = katss_init_kdata(kmer);
	for(uint64_t i=0; i<shuf->num_enrichments; i++) {
		double rval = prob->enrichments[i].enrichment / shuf->enrichments[i].enrichment;
		data->kmers[i].rval = opts->normalize ? log2(rval) : rval;
		data->kmers[i].kmer = (uint32_t)i;
	}

	katss_free_enrichments(prob);
exit:
	katss_free_enrichments(shuf);
	return data;

exit_error:
	katss_free_counter(test_counts);
	katss_free_counter(mono_counts);
	katss_free_counter(dint_counts);
	return NULL;
}

/**
 * @brief Compute the bootstrap enrichments of a dataset.
 * 
 * @param test Testfile to be used for computation
 * @param ctrl Control file to be used for computation
 * @param opts Options to modify output
 * @return KatssData* Data containing rval's, stdev, and pvalue
 */
static KatssData *
bootstrap_regular(const char *test, const char *ctrl, KatssOptions *opts)
{
	KatssCounter *test_counts = NULL;
	KatssCounter *ctrl_counts = NULL;
	KatssData *enrichments    = NULL;
	unsigned int seed         = opts->seed;
	unsigned int kmer         = opts->kmer;
	int sample                = opts->bootstrap_sample;
	int threads               = opts->threads;

	/* Create T-test aggregates */
	uint64_t total = 1ULL << (2*opts->kmer);
	t_test2_aggregate **ttest2 = s_malloc(sizeof *ttest2 * total);
	for(uint64_t i=0; i<total; i++)
		ttest2[i] = t_test2_create();

	/* Compute bootstrap values */
	double test_val, ctrl_val;
	for(int i=0; i<opts->bootstrap_iters; i++) {
		test_counts = katss_count_kmers_bootstrap_mt(test, kmer, sample, &seed, threads);
		ctrl_counts = katss_count_kmers_bootstrap_mt(ctrl, kmer, sample, &seed, threads);
		if(test_counts == NULL || ctrl_counts == NULL)
			goto exit_error;

		for(uint64_t k=0; k<total; k++) {
			katss_get_from_hash(test_counts, KATSS_DOUBLE, &test_val, (uint32_t)k);
			katss_get_from_hash(ctrl_counts, KATSS_DOUBLE, &ctrl_val, (uint32_t)k);
			test_val = test_val == 0 ? NAN : test_val;
			ctrl_val = ctrl_val == 0 ? NAN : ctrl_val;

			/* Update the t-test aggregate */
			t_test2_update(ttest2[k], test_val, ctrl_val);

			/* Use unused df and pval to be able to store rval stdev */
			if(!isnan(test_val) && !isnan(ctrl_val))
				running_stdev(test_val/ctrl_val, &ttest2[k]->df, &ttest2[k]->pval, i+1);
		}

		/* Free the counters */
		katss_free_counter(test_counts);
		katss_free_counter(ctrl_counts);
	}

	/* Finalize the bootstrap */
	enrichments = katss_init_kdata(opts->kmer);
	for(uint64_t i=0; i<enrichments->num_kmers; i++) {
		enrichments->kmers[i].kmer = i;
		enrichments->kmers[i].stdev = sqrt(ttest2[i]->pval / (opts->bootstrap_iters - 1));
		enrichments->kmers[i].rval = opts->normalize ? log2(ttest2[i]->df) : ttest2[i]->df;

		t_test2_finalize(ttest2[i]);
		enrichments->kmers[i].pval = ttest2[i]->pval;
		t_test2_destroy(ttest2[i]);
	}

	/* Clean up and return */
	free(ttest2);
	return enrichments;

exit_error:
	katss_free_counter(test_counts);
	katss_free_counter(ctrl_counts);
	for(uint64_t i=0; i<total; i++)
		t_test2_destroy(ttest2[i]);
	free(ttest2);
	return NULL;
}

/**
 * @brief Compute the bootstrap enrichments of using the probabilistic method.
 * 
 * @param test Testfile to be used for computation
 * @param opts Options to modify output
 * @return KatssData* Data containing rval's, stdev, and pvalue
 */
static KatssData *
bootstrap_probs(const char *test, KatssOptions *opts)
{
	KatssCounter *test_counts = NULL;
	KatssCounter *mono_counts = NULL;
	KatssCounter *dint_counts = NULL;
	KatssData *enrichments    = NULL;
	unsigned int kmer         = opts->kmer;
	int sample                = opts->bootstrap_sample;
	int threads               = opts->threads;
	unsigned int seed1,seed2,seed3;
	seed1 = seed2 = seed3 = opts->seed;

	/* Create T-test aggregates */
	uint64_t total = 1ULL << (2*opts->kmer);
	t_test2_aggregate **ttest2 = s_malloc(sizeof *ttest2 * total);
	for(uint64_t i=0; i<total; i++)
		ttest2[i] = t_test2_create();

	/* Compute bootstrap values */
	double test_val;
	for(int i=0; i<opts->bootstrap_iters; i++) {
		test_counts = katss_count_kmers_bootstrap_mt(test, kmer, sample, &seed1, threads);
		mono_counts = katss_count_kmers_bootstrap_mt(test, 1,    sample, &seed2, threads);
		dint_counts = katss_count_kmers_bootstrap_mt(test, 2,    sample, &seed3, threads);
		if(test_counts == NULL || mono_counts == NULL || dint_counts == NULL)
			goto exit_error;

		for(uint64_t k=0; k<total; k++) {
			katss_get_from_hash(test_counts, KATSS_DOUBLE, &test_val, (uint32_t)k);
			double ctrl_val = katss_predict_kmer_freq((uint32_t)k, kmer, mono_counts, dint_counts);
			ctrl_val *= katss_get_total(test_counts);

			/* Update the t-test aggregate */
			t_test2_update(ttest2[k], test_val, ctrl_val);

			/* Use unused df and pval to be able to store rval stdev */
			test_val /= katss_get_total(test_counts); // normalize the test_val
			ctrl_val = katss_predict_kmer_freq((uint32_t)k, kmer, mono_counts, dint_counts);
			running_stdev(test_val/ctrl_val, &ttest2[k]->df, &ttest2[k]->pval, i+1);
		}

		/* Free the counters */
		katss_free_counter(test_counts);
		katss_free_counter(mono_counts);
		katss_free_counter(dint_counts);
	}

	/* Finalize the bootstrap */
	enrichments = katss_init_kdata(opts->kmer);
	for(uint64_t i=0; i<enrichments->num_kmers; i++) {
		enrichments->kmers[i].kmer = i;
		enrichments->kmers[i].stdev = sqrt(ttest2[i]->pval / (opts->bootstrap_iters - 1));
		if(opts->normalize) {
			enrichments->kmers[i].rval = log2(ttest2[i]->df); // df holds rval
		} else {
			enrichments->kmers[i].rval = ttest2[i]->df; // df holds rval
		}

		t_test2_finalize(ttest2[i]);
		enrichments->kmers[i].pval = ttest2[i]->pval;
		t_test2_destroy(ttest2[i]);
	}

	/* Clean up and return */
	free(ttest2);
	return enrichments;

exit_error:
	katss_free_counter(test_counts);
	katss_free_counter(mono_counts);
	katss_free_counter(dint_counts);
	for(uint64_t i=0; i<total; i++)
		t_test2_destroy(ttest2[i]);
	free(ttest2);
	return NULL;
}

/**
 * @brief Compute the bootstrap enrichments of shuffled sequences.
 * 
 * @param test Testfile to be used for computation
 * @param opts Options to modify output
 * @return KatssData* Data containing rval's, stdev, and pvalue
 */
static KatssData *
bootstrap_ushuffle(const char *test, KatssOptions *opts)
{
	KatssCounter *test_counts = NULL;
	KatssCounter *shuf_counts = NULL;
	KatssData *enrichments    = NULL;
	unsigned int kmer         = opts->kmer;
	int klet                  = opts->probs_ntprec;
	int sample                = opts->bootstrap_sample;
	bool normalize            = opts->normalize;
	unsigned int seed1,seed2,seed3;
	seed1 = seed2 = seed3 = opts->seed;

	/* Create T-test aggregates */
	uint64_t total = 1ULL << (2*opts->kmer);
	t_test2_aggregate **ttest2 = s_malloc(sizeof *ttest2 * total);
	for(uint64_t i=0; i<total; i++)
		ttest2[i] = t_test2_create();

	/* Compute bootstrap values */
	for(int i=1; i<=opts->bootstrap_iters; i++) {
		test_counts = katss_count_kmers_bootstrap(test, kmer, sample, &seed1);
		if(test_counts == NULL)
			goto exit_error;
		shuf_counts = katss_count_kmers_ushuffle_bootstrap(test, kmer, klet, sample, &seed2);
		if(shuf_counts == NULL)
			goto exit_error;

		/* Update the statistics for all kmers in this iteration */
		for(uint64_t k=0; k<total; k++) {
			/* Obtain the predicted and actual counts for kmer k */
			double test_count, ctrl_count;
			katss_get_from_hash(test_counts, KATSS_DOUBLE, &test_count, (uint32_t)k);
			katss_get_from_hash(shuf_counts, KATSS_DOUBLE, &ctrl_count, (uint32_t)k);

			/* Update the t-test aggregate */
			t_test2_update(ttest2[k], test_count, ctrl_count);

			/* Use unused df and pval to be able to store rval stdev */
			test_count /= katss_get_total(test_counts);
			ctrl_count /= katss_get_total(shuf_counts);
			running_stdev(test_count/ctrl_count, &ttest2[k]->df, &ttest2[k]->pval, i);
		}

		/* Free the counters */
		katss_free_counter(test_counts);
		katss_free_counter(shuf_counts);
	}

	/* Finalize the bootstrap */
	enrichments = katss_init_kdata(opts->kmer);
	for(uint64_t i=0; i<enrichments->num_kmers; i++) {
		/* df currently holds the rval (due to running_stdev) */
		double rval = ttest2[i]->df;
		enrichments->kmers[i].rval = normalize ? log2(rval) : rval;
		enrichments->kmers[i].kmer = i;
		enrichments->kmers[i].stdev = sqrt(ttest2[i]->pval / (opts->bootstrap_iters - 1));

		/* Finalize the t-test to get pvalue */
		t_test2_finalize(ttest2[i]);
		enrichments->kmers[i].pval = ttest2[i]->pval;
		t_test2_destroy(ttest2[i]);
	}

	/* Clean up and return */
	free(ttest2);
	return enrichments;

exit_error:
	katss_free_counter(test_counts);
	katss_free_counter(shuf_counts);
	for(uint64_t i=0; i<total; i++)
		t_test2_destroy(ttest2[i]);
	free(ttest2);
	return NULL;
}


/**
 * @brief Compute the bootstraped enrichments using both shuffled and
 * probabilistic enrichments.
 * 
 * Enrichments are computed as the probabilistic enrichment of the shuffled
 * sequences over the probabilistic enrichment of the dataset. Each enrichment
 * is calculated using a subsample of the dataset (25% by default), over n
 * number of iterations. The standard deviation of the enrichments are
 * calculated using the enrichments (test_rval / shuffled_rval). Though, the
 * p-value is calculated though the two sample T-test, where sample 1 is all
 * the test_rval's, and sample 2 is shuffled_rval.
 * 
 * @param test Testfile to be used for computation
 * @param opts Options to modify output
 * @return KatssData* Data containing rval's, stdev, and pvalue
 */
static KatssData *
bootstrap_both(const char *test, KatssOptions *opts)
{
	KatssCounter *test_counts = NULL;
	KatssCounter *mono_counts = NULL;
	KatssCounter *dint_counts = NULL;
	KatssEnrichments *shuf    = NULL;
	KatssEnrichments *prob    = NULL;
	KatssData *enrichments    = NULL;
	unsigned int kmer         = opts->kmer;
	int klet                  = opts->probs_ntprec;
	int sample                = opts->bootstrap_sample;
	int threads               = opts->threads;
	unsigned int seed1,seed2,seed3,seed4,seed5,seed6;
	seed1 = seed2 = seed3 = seed4 = seed5 = seed6 = opts->seed;

	/* Create T-test aggregates */
	uint64_t total = 1ULL << (2*opts->kmer);
	t_test2_aggregate **ttest2 = s_malloc(sizeof *ttest2 * total);
	for(uint64_t i=0; i<total; i++)
		ttest2[i] = t_test2_create();

	/* Compute bootstrap values */
	for(int i=1; i<=opts->bootstrap_iters; i++) {
		/* Get the shuffled counts */
		test_counts = katss_count_kmers_ushuffle_bootstrap(test, kmer, klet, sample, &seed1);
		mono_counts = katss_count_kmers_ushuffle_bootstrap(test, 1,    klet, sample, &seed2);
		dint_counts = katss_count_kmers_ushuffle_bootstrap(test, 2,    klet, sample, &seed3);
		if(test_counts == NULL || mono_counts == NULL || dint_counts == NULL)
			goto exit_error;

		/* Compute probabilistic enrichmenst of shuffled counts */
		shuf = katss_compute_prob_enrichments(test_counts, mono_counts, dint_counts, false);
		if(shuf == NULL)
			goto exit_error;

		katss_free_counter(test_counts);
		katss_free_counter(mono_counts);
		katss_free_counter(dint_counts);

		/* Compute probabilistic enrichment of dataset */
		test_counts = katss_count_kmers_bootstrap_mt(test, kmer, sample, &seed4, threads);
		mono_counts = katss_count_kmers_bootstrap_mt(test, 1,    sample, &seed5, threads);
		dint_counts = katss_count_kmers_bootstrap_mt(test, 2,    sample, &seed6, threads);
		if(test_counts == NULL || mono_counts == NULL || dint_counts == NULL)
			goto exit_error_probs;
		prob = katss_compute_prob_enrichments(test_counts, mono_counts, dint_counts, false);
		if(prob == NULL)
			goto exit_error_probs;
		katss_free_counter(test_counts);
		katss_free_counter(mono_counts);
		katss_free_counter(dint_counts);

		/* Update the statistics for all kmers in this iteration */
		for(uint64_t k=0; k<total; k++) {
			/* Update the t-test aggregate */
			double test_rval = prob->enrichments[k].enrichment;
			double ctrl_rval = shuf->enrichments[k].enrichment;
			t_test2_update(ttest2[k], test_rval, ctrl_rval);

			/* Store the standard deviation for R value */
			running_stdev(test_rval/ctrl_rval, &ttest2[k]->df, &ttest2[k]->pval, i);
		}

		/* Free the enrichments */
		katss_free_enrichments(prob);
		katss_free_enrichments(shuf);
	}

	/* Finalize the bootstrap */
	enrichments = katss_init_kdata(opts->kmer);
	for(uint64_t i=0; i<enrichments->num_kmers; i++) {
		enrichments->kmers[i].kmer = i;
		enrichments->kmers[i].stdev = sqrt(ttest2[i]->pval / (opts->bootstrap_iters - 1));
		if(opts->normalize) {
			enrichments->kmers[i].rval = log2(ttest2[i]->df);
		} else {
			enrichments->kmers[i].rval = ttest2[i]->df; // df holds rval currently
		} 

		t_test2_finalize(ttest2[i]);
		enrichments->kmers[i].pval = ttest2[i]->pval;
		t_test2_destroy(ttest2[i]);
	}

	/* Clean up and return */
	free(ttest2);
	return enrichments;

exit_error_probs:
	katss_free_enrichments(shuf);
exit_error:
	katss_free_counter(test_counts);
	katss_free_counter(mono_counts);
	katss_free_counter(dint_counts);
	for(uint64_t i=0; i<total; i++)
		t_test2_destroy(ttest2[i]);
	free(ttest2);
	return NULL;
}

KatssData *
katss_enrichment(const char *test, const char *ctrl, KatssOptions *opts)
{
	/* Make sure test file was passed */
	if(test == NULL)
		return NULL;

	/* Parse the options */
	if(katss_parse_options(opts) != 0)
		return NULL;

	/* Ensure control file exists if not using a probabilistic algo */
	if(ctrl == NULL && opts->probs_algo == KATSS_PROBS_NONE && opts->enable_warnings)
		error_message("katss_enrichment: If no probabilistic algorithm is set,"
		              " `ctrl' can't be NULL");
	if(ctrl == NULL && opts->probs_algo == KATSS_PROBS_NONE)
		return NULL;
	if(ctrl && opts->probs_algo != KATSS_PROBS_NONE && opts->enable_warnings)
		warning_message("katss_enrichment: Ignoring `ctrl=(%s)'",ctrl);

	/* BEGIN COMPUTATION: No bootstrap */
	KatssData *data = NULL;
	if(opts->bootstrap_iters == 0) {
		switch(opts->probs_algo) {
		case KATSS_PROBS_NONE:     data = regular(test, ctrl, opts); break;
		case KATSS_PROBS_REGULAR:  data = probs(test, opts);         break;
		case KATSS_PROBS_USHUFFLE: data = ushuffle(test, opts);      break;
		case KATSS_PROBS_BOTH:     data = both(test, opts);          break;
		default: return NULL; // silence compiler warnings
		}

	/* BEGIN COMPUTATION: bootstrap */
	} else {
		switch(opts->probs_algo) {
		case KATSS_PROBS_NONE:     data = bootstrap_regular(test, ctrl, opts); break;
		case KATSS_PROBS_REGULAR:  data = bootstrap_probs(test, opts);         break;
		case KATSS_PROBS_USHUFFLE: data = bootstrap_ushuffle(test, opts);      break;
		case KATSS_PROBS_BOTH:     data = bootstrap_both(test, opts);          break;
		default: return NULL; // silence compiler warnings
		}
	}

	/* If data is NULL, return NULL */
	if(data == NULL)
		return NULL;

	/* Sort if necessary */
	if(opts->sort_enrichments)
		qsort(data->kmers, data->num_kmers, sizeof *data->kmers, compare);

	/* Return data! */
	return data;
}
