#ifndef KATSS_ENRICHMENTS_H
#define KATSS_ENRICHMENTS_H

#include <stdint.h>
#include <stdbool.h>
#include "counter.h"

#ifdef __cplusplus
extern "C" {
#endif 

typedef struct KatssEnrichment {
	double enrichment;
	uint32_t key;
} KatssEnrichment;

typedef struct KatssEnrichments {
	KatssEnrichment *enrichments;
	uint64_t num_enrichments;
} KatssEnrichments;

/* Enrichments functions */
KatssEnrichments *katss_compute_enrichments(KatssCounter *test, KatssCounter *control, bool normalize);
KatssEnrichments *katss_enrichments(const char *test_file, const char *control_file, unsigned int kmer, bool normalize);
KatssEnrichments *katss_compute_prob_enrichments(KatssCounter *test, KatssCounter *mono, 
                                                 KatssCounter *dint, bool normalize);
KatssEnrichments *katss_prob_enrichments(const char *test_file, unsigned int kmer, bool normalize);

/* IKKE Functions */
KatssEnrichments *katss_ikke_mt(const char *test_file, const char *control_file, unsigned int kmer, 
                                uint64_t iterations, bool normalize, int threads);
KatssEnrichments *katss_ikke_(const char *test_file, const char *control_file, unsigned int kmer, uint64_t iterations, bool normalize);
KatssEnrichments *katss_prob_ikke_mt(const char *test_file, unsigned int kmer, uint64_t iterations, bool normalize, int threads);
KatssEnrichments *katss_prob_ikke(const char *test_file, unsigned int kmer, uint64_t iterations, bool normalize);

KatssEnrichment katss_top_enrichment(KatssCounter *test, KatssCounter *control, bool normalize);
KatssEnrichment katss_top_prediction(KatssCounter *test, KatssCounter *mono, KatssCounter *dint, bool normalize);
void katss_free_enrichments(KatssEnrichments *enrichments);
void katss_sort_enrichments(KatssEnrichments *enrichments);

#ifdef __cplusplus
}
#endif

#endif // KATSS_ENRICHMENTS_H
