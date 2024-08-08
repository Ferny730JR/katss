#include <string.h>
#include <math.h>

#include <ViennaRNA/fold.h>
#include <ViennaRNA/part_func.h>

#include "rnafiles.h"
#include "memory_utils.h"
#include "bpp_tables.h"
#include "structure.h"
#include "string_utils.h"

#if (__STDC_NO_THREADS__)
#  include "tinycthread.h"
#else
#  include <threads.h>
#endif

#ifdef BUFFER_SIZE
#undef BUFFER_SIZE
#endif

#define BUFFER_SIZE 65536

struct record_data {
	char *sequence;
	kmerHashTable *counts_table;
	BppOptions *opts;
	RnaFile read_file;
};

typedef struct record_data record_data;

/* Helper functions */
static bool
is_nucleotide(char character)
{
	switch(character) {
		case 'A':   return true;
		case 'a':   return true;
		case 'C':   return true;
		case 'c':   return true;
		case 'G':   return true;
		case 'g':   return true;
		case 'T':   return true;
		case 't':   return true;
		case 'U':   return true;
		case 'u':   return true;
		default:    return false;
	}
}

static char
determine_filetype(const char *file)
{
	/* Open the RnaFile, return 'e' upon error */
	RnaFile reads_file = rnafopen(file, "b");
	if(reads_file == NULL) {
		char errbuf[1000];
		rnafstrerror_r(rnaferrno, errbuf, sizeof(errbuf));
		error_message("katss: %s: %s\n", file, errbuf);
		rnafclose(reads_file);
		return 'N';
	}

	char buffer[BUFFER_SIZE];
	int lines_read = 0;
	int fastq_score_lines = 0;
	int fasta_score_lines = 0;
	int sequence_lines = 0;

	while (rnafgets(reads_file, buffer, BUFFER_SIZE) != NULL && lines_read < 10) {
		lines_read++;
		char first_char = buffer[0];

		/* Check if the first line starts with '@' for FASTQ */
		if (first_char == '@' && lines_read % 4 == 1) {
			fastq_score_lines++;

		/* Check if the third line starts with '+' for FASTQ */
		} else if (first_char == '+' && lines_read % 4 == 3) {
			fastq_score_lines++;

		/* Check if the line starts with '>' or ';' for FASTA */
		/* TODO: fastq can potentially have a '>' or ';' in its quality score,
		meaning that file can be incorrectly guessed as fasta when it is fastq */
		} else if (first_char == '>' || first_char == ';') {
			fasta_score_lines++;
		} else {
			// Check for nucleotide characters
			int num_total = 0, num = 0;
			for(int i = 0; buffer[i] != '\0'; i++) {
				if(is_nucleotide(buffer[i])) {
					num++;
				}
				num_total++;
			}
			if((double)num/num_total > 0.9) {
				sequence_lines++;
			}
		}
	}
    rnafclose(reads_file);

    if (fastq_score_lines >= 2) {
        return 'q'; // fastq file
	} else if (fasta_score_lines > 0) {
		return 'a';
    } else if (sequence_lines == 10) {
        return 's'; // raw sequences file
    } else {
		error_message("Unable to read sequence from file.\nCurrent supported file types are:"
		              " FASTA, FASTQ, and file containing sequences per line.");
        return 'e'; // unsupported file type
    }
}

static RnaFile
rnafopen_detect(const char *filename)
{
	char filetype = determine_filetype(filename);
	if(filetype == 'N' || filetype == 'e')
		return NULL;
	char mode[2] = {filetype, '\0'};
	return rnafopen(filename, mode);
}

static float *
getPositionalProbabilities(char *sequence)
{
	vrna_ep_t *ptr, *pair_probabilities = NULL;
	float *positional_probabilities = s_calloc(strlen(sequence), sizeof(float));
	float  probability;

	/* Get the pair probabilities */
	vrna_pf_fold(sequence, NULL, &pair_probabilities);

	/* Move pair probabilities into array */
	for(ptr = pair_probabilities; ptr->i != 0; ptr++) {
		probability = ptr->p;
		positional_probabilities[ptr->i-1]+=probability;
		positional_probabilities[ptr->j-1]+=probability;
	}

	/* Clean up memory */
	free(pair_probabilities);

	return positional_probabilities;
}

static void
process_sequence(record_data *record)
{
	BppOptions *opts = record->opts;

	char  *sequence  = record->sequence;
	int   num_kmers_in_seq = strlen(sequence) - opts->kmer + 1;

	float *positional_probabilities = getPositionalProbabilities(sequence);

	/* Count kmers and their associated base-pair probability */
	for(int i=0; i<num_kmers_in_seq; i++) {
		int shift = opts->kmer+i;
		char tmp = *(sequence+shift);
		*(sequence+shift) = '\0'; // terminate the string to k-mer length

		kmer_add_value(record->counts_table, sequence+i, 1, opts->kmer);		
		/* Loop through bpp values in file */
		for(int j=i; j<opts->kmer+i; j++) {
			kmer_add_value(record->counts_table, sequence+i, positional_probabilities[j], j-i);
		}
		*(sequence+shift) = tmp;
	}

	free(positional_probabilities);
}

void
process_windows(record_data *record)
{
	char       *sequence  = record->sequence;
	BppOptions *opts       = record->opts;
	kmerHashTable *counts_table = record->counts_table;

	char    *window_seq;
	float   *window_probabilities;
	float   mean_probability;
	float   sum_probability;
	int     num_windows;
	int     seq_length;
	int     count_probs;

	seq_length = strlen(sequence);
	num_windows = seq_length - opts->window_size + 1;
	if(num_windows < 1) {
		process_sequence(record);
		return;
	}

	/* Initialize probability matrix with -1 */
	float probability_matrix[num_windows][seq_length];
	for(int row=0; row<num_windows; row++){
		for(int col=0; col<seq_length; col++) {
			probability_matrix[row][col]=-1;
		}
	}

	/* Fill the probability matrix with probabilities */
	for(int i = 0; i<num_windows; i++) {
		/* Get the window sequence */
		window_seq = sequence + i;
		char tmp = window_seq[opts->window_size];
		window_seq[opts->window_size] = '\0';

		/* Get probabilities from the window sequence */
		window_probabilities = getPositionalProbabilities(window_seq);
		window_seq[opts->window_size] = tmp;

		/* Dump probabilities to matrix*/
		for(int j = 0; j < opts->window_size; j++) {
			probability_matrix[i][j+i] = window_probabilities[j];
		}

		free(window_probabilities);
	}

	// Get the average of each column in probability matrix
	float positional_probabilities[seq_length];
	for(int col=0; col<seq_length; col++) {
		sum_probability = 0;
		count_probs     = 0;

		for(int row=0; row<num_windows; row++) {
			if(probability_matrix[row][col] == -1) {
				continue;
			}
			sum_probability+=probability_matrix[row][col];
			count_probs++;
		}

		mean_probability = sum_probability/count_probs;
		positional_probabilities[col] = mean_probability;
	}

	/* Fill counts_table with positional probabilities */
	int num_kmers_in_seq = seq_length - opts->kmer + 1;
	for(int i=0; i<num_kmers_in_seq; i++) {
		int shift = opts->kmer + i;
		char tmp = *(sequence + shift);
		*(sequence + shift) = '\0';
		kmer_add_value(counts_table, sequence+i, 1, opts->kmer);
		for(int j=i; j<opts->kmer+i; j++) {
			kmer_add_value(counts_table, sequence+i, positional_probabilities[j], j-i);
		}
		*(sequence + shift) = tmp;
	}
}

static void
process_record(record_data *record)
{
	/* If --seq-windows was provided, use the sliding window algorithm for calculations */
	if(record->opts->seq_windows) {
		process_windows(record);

	/* Default algorithm for getting BPP frequencies */
	} else {
		process_sequence(record);
	}
}

static int
bpp_kmer_count(void *arg)
{
	record_data *record = (record_data *)arg;
	char *sequence = s_malloc(BUFFER_SIZE * sizeof *sequence);

	record->sequence = sequence;
	while(true) {
		if(rnafgets(record->read_file, sequence, BUFFER_SIZE) == NULL)
			break;
		clean_seq(sequence, true);
		process_record(record);
	}
	free(sequence);
	return 0;
}

static kmerHashTable *
bpp_kmer_frequency(const char *filename, BppOptions *opts)
{
	kmerHashTable   *counts_table;
	RnaFile         read_file;

	read_file    = rnafopen_detect(filename);
	if(read_file == NULL)
		return NULL;
	counts_table = init_bpp_table(opts->kmer);

	/* Multi-threaded bpp counting */
	if(opts->threads > 1) {
		record_data *rd = s_malloc(opts->threads * sizeof *rd);
		thrd_t *jobs = s_malloc(opts->threads * sizeof *jobs);
		for(int i=0; i<opts->threads; i++) {
			rd[i].read_file = read_file;
			rd[i].opts = opts;
			rd[i].counts_table = counts_table;
			rd[i].sequence = NULL;
			thrd_create(&jobs[i], bpp_kmer_count, &rd[i]);
		}
		for(int i=0; i<opts->threads; i++) {
			thrd_join(jobs[i], NULL);
		}
		free(jobs);
		free(rd);
	/* Single-threaded bpp counting */
	} else {
		record_data *record = s_malloc(sizeof *record);
		record->sequence = NULL;
		record->counts_table = counts_table;
		record->opts = opts;
		record->read_file = read_file;
		bpp_kmer_count((void *)record);
		free(record);
	}

	rnafclose(read_file);

	/* Calculate the frequencies */
	int num_columns = counts_table->cols-1;
	for(int i=0; i<counts_table->capacity; i++) {
		if( !counts_table->entries[i] || counts_table->entries[i]->values[num_columns]<1) {
			continue;
		}
		for(int j=0; j<num_columns; j++) {
			int total_count=(int)counts_table->entries[i]->values[num_columns];
			counts_table->entries[i]->values[j]/=total_count;
		}
	}

	return counts_table;
}

static int
bpp_compare(const void *p1, const void *p2)
{
	const Entry *entryA = *(Entry **)p1;
	const Entry *entryB = *(Entry **)p2;

	if(!entryA && !entryB)
		return 0;
	if(!entryA)
		return 1;
	if(!entryB)
		return -1;

	int mean_index = entryA->num_values-1;
	if(entryA->values[mean_index] > entryB->values[mean_index])
		return -1;
	if(entryA->values[mean_index] < entryB->values[mean_index])
		return 1;
	return 0;
}

static kmerHashTable *
bpp_enrichment(kmerHashTable *control_frq, kmerHashTable *bound_frq, int kmer)
{
	kmerHashTable   *enrichments_table;
	double          enrichment;
	double          *bound_values;
	double          *control_values;
	double          *enrichment_values;
	double          mean_enrichment;

	enrichments_table = init_bpp_table(kmer);

	// Get the log2 fold change for each kmer
	for(size_t i = 0; i < bound_frq->capacity; i++) {
		if(bound_frq->entries[i] == NULL || control_frq->entries[i] == NULL) {
			continue;
		}

		bound_values   = bound_frq->entries[i]->values;
		control_values = control_frq->entries[i]->values;

		enrichments_table->entries[i] = s_malloc(sizeof(Entry));
		enrichments_table->entries[i]->num_values = enrichments_table->cols;
		enrichments_table->entries[i]->hash = i;
		enrichments_table->entries[i]->values = s_calloc(enrichments_table->cols, sizeof(double));

		for(int j = 0; j < kmer; j++) {
			if(bound_values[j] == 0. || control_values[j] == 0.) {
				enrichment = 0.;
			} else {
				enrichment = log2(bound_values[j]/control_values[j]);
			}
			// todo: create kmer_set_value function
			enrichments_table->entries[i]->values[j] = enrichment;
		}
	}

	for(size_t i = 0; i < enrichments_table->capacity; i++) {
		if(enrichments_table->entries[i] == NULL) {
			continue;
		}

		enrichment_values = enrichments_table->entries[i]->values;

		mean_enrichment = 0.0f;
		for(int j = 0; j < kmer; j++) {
			mean_enrichment += enrichment_values[j];
		}

		mean_enrichment/=kmer;
		enrichments_table->entries[i]->values[kmer] = mean_enrichment;
	}

	qsort(enrichments_table->entries, enrichments_table->capacity,
			sizeof(*enrichments_table->entries), bpp_compare);

	return enrichments_table;
}

void
katss_free_bpp(kmerHashTable *bpp)
{
	free_kmer_table(bpp);
}

void
katss_bpp_init_default_opts(BppOptions *opts)
{
	opts->kmer = 5;
	opts->seq_windows = false;
	opts->threads = 1;
}

kmerHashTable *
katss_bpp(const char *test_file, const char *ctrl_file, BppOptions *opts)
{
	bool provided_opts = opts != NULL;
	if(opts == NULL) {
		opts = s_malloc(sizeof *opts);
		katss_bpp_init_default_opts(opts);
	}

	kmerHashTable *enrichments = NULL;

	if(test_file == NULL || ctrl_file == NULL)
		return NULL;
	if(opts->kmer > 16)
		return NULL;
	opts->threads = MAX2(opts->threads, 1);
	opts->threads = MIN2(opts->threads, 128);

	kmerHashTable *test_frq = bpp_kmer_frequency(test_file, opts);
	if(test_frq == NULL)
		goto exit;
	kmerHashTable *ctrl_frq = bpp_kmer_frequency(ctrl_file, opts);
	if(ctrl_frq == NULL)
		goto undo_test;
	enrichments = bpp_enrichment(ctrl_frq, test_frq, opts->kmer);

	if(!provided_opts)
		free(opts);
	free_kmer_table(ctrl_frq);
undo_test:
	free_kmer_table(test_frq);
exit:
	return enrichments;
}
