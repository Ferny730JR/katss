#include <stddef.h>
#include <stdlib.h>

#include "katss.h"
#include "katss_helpers.h"
#include "memory_utils.h"

#include "enrichments.h"

static KatssData *
regular(const char *test, const char *ctrl, KatssOptions *opts)
{
	/* Compute iterative kmer knockout enrichments */
	KatssEnrichments *enr;
	enr = katss_ikke_mt(test, ctrl, opts->kmer, opts->iters, opts->normalize, opts->threads);
	if(enr == NULL)
		return NULL;

	/* Move enrichments to data */
	KatssData *data;
	if((data = katss_init_kdata(opts->kmer)) == NULL)
		goto exit;
	for(uint64_t i=0; i<enr->num_enrichments; i++) {
		data->kmers[i].kmer = enr->enrichments[i].key;
		data->kmers[i].rval = enr->enrichments[i].enrichment;
	}

exit:
	katss_free_enrichments(enr);
	return data;
}

static KatssData *
probs(const char *test, KatssOptions *opts)
{
	KatssEnrichments *enr;
	enr = katss_prob_ikke_mt(test, opts->kmer, opts->iters, opts->normalize, opts->threads);
	if(enr == NULL)
		return NULL;
	
	KatssData *data;
	if((data = katss_init_kdata(opts->kmer)) == NULL)
		goto exit;
	for(uint64_t i=0; i<enr->num_enrichments; i++) {
		data->kmers[i].kmer = enr->enrichments[i].key;
		data->kmers[i].rval = enr->enrichments[i].enrichment;
	}

exit:
	katss_free_enrichments(enr);
	return data;
}

static KatssData *
ushuffle(const char *test, KatssOptions *opts)
{
	KatssEnrichments *enr;
	enr = katss_ikke_shuffle(test, NULL, opts->kmer, opts->probs_ntprec, opts->iters, opts->normalize);
	if(enr == NULL)
		return NULL;
	
	KatssData *data = katss_init_kdata(opts->kmer);
	if(data == NULL)
		goto exit;
	for(uint64_t i=0; i<enr->num_enrichments; i++) {
		data->kmers[i].kmer = enr->enrichments[i].key;
		data->kmers[i].rval = enr->enrichments[i].enrichment;
	}

exit:
	katss_free_enrichments(enr);
	return data;
}

static KatssData *
both(const char *test, KatssOptions *opts)
{
	return NULL;
}

static KatssData *
bootstrap_regular(const char *test, const char *ctrl, KatssOptions *opts)
{
	return NULL;
}

static KatssData *
bootstrap_probs(const char *test, KatssOptions *opts)
{
	return NULL;
}

static KatssData *
bootstrap_ushuffle(const char *test, KatssOptions *opts)
{
	return NULL;
}

static KatssData *
bootstrap_both(const char *test, KatssOptions *opts)
{
	return NULL;
}

KatssData *
katss_ikke(const char *test, const char *ctrl, KatssOptions *opts)
{
	/* Make sure test file was passed */
	if(test == NULL)
		return NULL;

	/* Parse the options */
	if(katss_parse_options(opts) != 0)
		return NULL;
	
	/* Ensure control file exists if not using a probabilistic algo */
	if(ctrl == NULL && opts->probs_algo == KATSS_PROBS_NONE && opts->enable_warnings)
		error_message("katss_ikke: If no probabilistic algorithm is set,"
		              " `ctrl' can't be NULL");
	if(ctrl == NULL && opts->probs_algo == KATSS_PROBS_NONE)
		return NULL;
	if(ctrl && opts->probs_algo != KATSS_PROBS_NONE && opts->enable_warnings)
		warning_message("katss_enrichment: Ignoring `ctrl=(%s)'",ctrl);

	/* BEGIN COMPUTATION: No bootstrap */
	if(opts->bootstrap_iters == 0) {
		switch(opts->probs_algo) {
		case KATSS_PROBS_NONE:     return regular(test, ctrl, opts);
		case KATSS_PROBS_REGULAR:  return probs(test, opts);
		case KATSS_PROBS_USHUFFLE: return ushuffle(test, opts);
		case KATSS_PROBS_BOTH:     return both(test, opts);
		default: return NULL;
		}

	/* BEGIN COMPUTATION: bootstrap */
	} else {
		switch(opts->probs_algo) {
		case KATSS_PROBS_NONE:     return bootstrap_regular(test, ctrl, opts);
		case KATSS_PROBS_REGULAR:  return bootstrap_probs(test, opts);
		case KATSS_PROBS_USHUFFLE: return bootstrap_ushuffle(test, opts);
		case KATSS_PROBS_BOTH:     return bootstrap_both(test, opts);
		default: return NULL;
		}
	}
}
