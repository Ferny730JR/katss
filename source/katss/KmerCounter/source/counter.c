#include <stdbool.h>
#include <errno.h>
#include <string.h>

#if (__STDC_NO_THREADS__)
#  include "tinycthread.h"
#else
#  include <threads.h>
#endif

#include "counter.h"
#include "hash_functions.h"
#include "memory_utils.h"
#include "rnafiles.h"

#define BUFFER_SIZE 65536U

struct threadinfo {
	RnaFile rnafile;
	KatssCounter *counter;
	unsigned int kmer;
	char filetype;
};
typedef struct threadinfo threadinfo;

/*============ Counting Function Declarations ============*/
static KatssCounter *
count_file(const char *filename, unsigned int kmer, const char filetype);
static int
count_file_mt(void *arg);

/*============= Helper Function Declarations =============*/
static char
determine_filetype(const char *filename);

static bool
is_nucleotide(char character);

/*============= Actual Functions Declarations =============*/
KatssCounter *
katss_count_kmers(const char *filename, unsigned int kmer)
{
	char filetype = determine_filetype(filename);
	if(filetype == 'e') { /* Error determining filetype */
		return NULL;
	} else if (filetype == 'N') { /* Error opening file */
		return NULL;
	}

	KatssCounter *counter = count_file(filename, kmer, filetype);
	return counter;
}


KatssCounter *
katss_count_kmers_mt(const char *filename, unsigned int kmer, int threads)
{
	char filetype = determine_filetype(filename);
	if(filetype == 'e') { /* Error determining filetype */
		return NULL;
	} else if (filetype == 'N') { /* Error opening file */
		return NULL;
	}

	threads = MAX2(threads, 1);
	threads = MIN2(threads, 128);

	/* Open RnaFile for reading */
	char mode[2] = { 0 };
	mode[0] = filetype == 'r' ? 's' : filetype;
	RnaFile file = rnafopen(filename, mode);
	if(file == NULL) {
		char *error = rnafstrerror(rnaferrno);
		warning_message("rnafopen: error %d: %s",rnaferrno,error);
		free(error);
		return NULL;
	}

	/* Initialize counter */
	KatssCounter *counter = katss_init_counter(kmer);
	if(counter == NULL) {
		rnafclose(file);
		return NULL;
	}

	threadinfo *jobarg = s_malloc(threads * sizeof *jobarg);
	thrd_t *jobs = s_malloc(threads * sizeof *jobs);
	for(int i=0; i<threads; i++) {
		jobarg[i].rnafile = file;
		jobarg[i].counter = counter;
		jobarg[i].kmer = kmer;
		jobarg[i].filetype = filetype;

		/* Start threads */
		thrd_create(&jobs[i], count_file_mt, &jobarg[i]);
	}

	for(int i=0; i<threads; i++) {
		thrd_join(jobs[i], NULL);
	}

	/* Free resources */
	rnafclose(file);
	free(jobs);
	free(jobarg);

	return counter;
}


static KatssCounter *
count_file(const char *filename, unsigned int kmer, const char filetype)
{
	KatssCounter *counter = NULL;

	/* Open file and prepare counter & hasher */
	RnaFile read_file = rnafopen(filename, "b");
	if(read_file == NULL)
		goto exit;

	KatssHasher *hasher = katss_init_hasher(kmer, filetype);
	if(hasher == NULL)
		goto cleanup_file;

	counter = katss_init_counter(kmer);
	if(counter == NULL)
		goto cleanup_hasher;

	/* Prepare file reading & hash int */
	char buffer[BUFFER_SIZE+1] = { 0 };
	size_t still_reading;
	uint32_t hash_value;

	do {
		still_reading = rnafread_unlocked(read_file, buffer, BUFFER_SIZE);
		buffer[still_reading] = '\0';

		katss_set_seq(hasher, buffer, filetype);
		while(katss_get_fh(hasher, &hash_value, filetype)) {
			katss_increment(counter, hash_value);
		}
	} while(still_reading == BUFFER_SIZE);

	/* If error was encountered while reading report and return NULL */
	if(still_reading == 0 && rnaferrno) {
		katss_free_counter(counter);
		counter = NULL;
		char *err = rnafstrerror(rnaferrno);
		error_message("katss: %d: %s", rnaferrno, err);
		free(err);
	}

cleanup_hasher:
	free(hasher);
cleanup_file:
	rnafclose(read_file);
exit:
	return counter;
}

static int
count_file_mt(void *arg)
{
	threadinfo *args = (threadinfo *)arg;
	char *buffer = s_malloc(BUFFER_SIZE * sizeof *buffer);

	KatssHasher *hasher = katss_init_hasher(args->kmer, '\0');
	if(hasher == NULL) {
		return 1;
	}

	/* Megabyte to store counts */
	size_t num_counts = 250000;
	uint32_t *hash_values = s_malloc(num_counts * sizeof *hash_values);
	size_t cur_hash = 0;

	/* Begin counting */
	while(rnafread(args->rnafile, buffer, BUFFER_SIZE)) {
		katss_set_seq(hasher, buffer, args->filetype);
		while(katss_get_fh(hasher, &hash_values[cur_hash], args->filetype)) {
			if(++cur_hash == num_counts) { // begin flushing
				katss_increments(args->counter, hash_values, cur_hash);
				cur_hash = 0;
			}
		}
	}

	/* Flush values */
	katss_increments(args->counter, hash_values, cur_hash);

	/* Free resources */
	free(hasher);
	free(buffer);
	free(hash_values);

	return 0;
}
/*==============================================================
|  Helper Functions                                            |
==============================================================*/
static char
determine_filetype(const char *file)
{
	/* Open the RnaFile, return 'e' upon error */
	RnaFile reads_file = rnafopen(file, "b");
	if(reads_file == NULL) {
		error_message("katss: %s: %s", file, strerror(errno));
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
        return 'r'; // raw sequences file
    } else {
		error_message("Unable to read sequence from file.\nCurrent supported file types are:"
		              " FASTA, FASTQ, and file containing sequences per line.");
        return 'e'; // unsupported file type
    }
}


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
