#ifdef _WIN32
#define _CRT_RAND_S
#define rand_r rand_s
#endif
#include <stdlib.h> // for random functions 
#include <stdbool.h>
#include <errno.h>
#include <string.h>
#include <time.h>

#if (__STDC_NO_THREADS__)
#  include "tinycthread.h"
#else
#  include <threads.h>
#endif

#include "counter.h"
#include "hash_functions.h"
#include "memory_utils.h"
#include "ushuffle.h"
#include "seqfile.h"
#include "thread_safe_rand.h"
#define BUFFER_SIZE 65536U

struct threadinfo {
	SeqFile seqfile;
	KatssCounter *counter;
	unsigned int kmer;
	int sample;
	unsigned int *seed;
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
	/* Threads should be at least one, and at most 128 */
	threads = MAX2(threads, 1);
	threads = MIN2(threads, 128);

	/* If one thread, use single threaded computation */
	if(threads == 1)
		return katss_count_kmers(filename, kmer);

	/* Begin multithreaded computations */
	char filetype = determine_filetype(filename);
	if(filetype == 'e') { /* Error determining filetype */
		return NULL;
	} else if (filetype == 'N') { /* Error opening file */
		return NULL;
	}

	/* Open SeqFile for reading */
	char mode[2] = { 0 };
	mode[0] = filetype == 'r' ? 's' : filetype;
	SeqFile file = seqfopen(filename, mode);
	if(file == NULL) {
		warning_message("seqfopen: error %d: %s",seqferrno,seqfstrerror(seqferrno));
		return NULL;
	}

	/* Initialize counter */
	KatssCounter *counter = katss_init_counter(kmer);
	if(counter == NULL) {
		seqfclose(file);
		return NULL;
	}

	threadinfo *jobarg = s_malloc(threads * sizeof *jobarg);
	thrd_t *jobs = s_malloc(threads * sizeof *jobs);
	for(int i=0; i<threads; i++) {
		jobarg[i].seqfile = file;
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
	seqfclose(file);
	free(jobs);
	free(jobarg);

	return counter;
}


KatssCounter *
katss_count_kmers_bootstrap(const char *filename, unsigned int kmer,
                            int sample, unsigned int *seed)
{
	char filetype = determine_filetype(filename);
	if(filetype == 'e' || filetype == 'N')
		return NULL;

	KatssCounter *counter = NULL;

	/* Initialize buffer */
	char *buffer = s_calloc(BUFFER_SIZE, sizeof *buffer);
	uint32_t hash_value;

	/* Open file and hasher */
	buffer[0] = filetype == 'r' ? 'b' : filetype;
	SeqFile read_file = seqfopen(filename, buffer);
	if(read_file == NULL)
		goto exit;
	
	KatssHasher *hasher = katss_init_hasher(kmer, filetype);
	if(hasher == NULL)
		goto cleanup_file;
	
	counter = katss_init_counter(kmer);
	if(counter == NULL)
		goto cleanup_hasher;

	/* sample should be between 1-100000 */
	sample = MAX2(sample, 1);
	sample = MIN2(sample, 100000);

	unsigned int local_seed;
	if(seed == NULL) {
		local_seed = time(NULL);
		seed = &local_seed;
	}

	while(seqfgets_unlocked(read_file, buffer, BUFFER_SIZE)) {
		if(rand_r(seed) % 100000 >= sample)
			continue;
		katss_set_seq(hasher, buffer, filetype);
		while(katss_get_fh(hasher, &hash_value, filetype)) {
			katss_increment(counter, hash_value);
		}
	}

	if(seqferrno) {
		error_message("katss: sample: %s\n", seqfstrerror_r(seqferrno, buffer, BUFFER_SIZE));
		katss_free_counter(counter);
		counter = NULL;
	}

	free(buffer);
cleanup_hasher:
	free(hasher);
cleanup_file:
	seqfclose(read_file);
exit:
	return counter;
}


static int
count_file_bootstrap_mt(void *arg)
{
	threadinfo *args = (threadinfo *)arg;
	char *buffer = s_malloc(BUFFER_SIZE * sizeof *buffer);
	thread_safe_rand_t *tsr = thread_safe_rand_init();

	KatssHasher *hasher = katss_init_hasher(args->kmer, '\0');
	if(hasher == NULL)
		return 1;

	/* Megabyte to store counts */
	size_t num_counts = 250000;
	uint32_t *hash_values = s_malloc(num_counts * sizeof *hash_values);
	size_t cur_hash = 0;

	/* Begin counting */
	while(seqfgets(args->seqfile, buffer, BUFFER_SIZE)) {
		if(thread_safe_rand_r(tsr, args->seed) % 100000 >= args->sample)
			continue;
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
	thread_safe_rand_free(tsr);
	free(hasher);
	free(buffer);
	free(hash_values);

	return 0;
}


KatssCounter *
katss_count_kmers_bootstrap_mt(const char *filename, unsigned int kmer,
                               int sample, unsigned int *seed, int threads)
{
	threads = MAX2(threads, 1);
	threads = MIN2(threads, 128);

	/* Process single-threaded computation */
	if(threads == 1)
		return katss_count_kmers_bootstrap(filename, kmer, sample, seed);

	/* Multi-threading */
	char filetype = determine_filetype(filename);
	if(filetype == 'e' || filetype == 'N')
		return NULL;

	/* sample should be between 1-100000 */
	sample = MAX2(sample, 1);
	sample = MIN2(sample, 100000);

	/* Open SeqFile for reading */
	char mode[2] = { 0 };
	mode[0] = filetype == 'r' ? 's' : filetype;
	SeqFile file = seqfopen(filename, mode);
	if(file == NULL) {
		warning_message("seqfopen: error %d: %s",seqferrno,seqfstrerror(seqferrno));
		return NULL;
	}

	/* Initialize counter */
	KatssCounter *counter = katss_init_counter(kmer);
	if(counter == NULL) {
		seqfclose(file);
		return NULL;
	}

	threadinfo *jobarg = s_malloc(threads * sizeof *jobarg);
	thrd_t *jobs = s_malloc(threads * sizeof *jobs);
	for(int i=0; i<threads; i++) {
		jobarg[i].seqfile = file;
		jobarg[i].counter = counter;
		jobarg[i].kmer = kmer;
		jobarg[i].filetype = filetype;
		jobarg[i].sample = sample;
		jobarg[i].seed = seed;

		/* Start threads */
		thrd_create(&jobs[i], count_file_bootstrap_mt, &jobarg[i]);
	}

	for(int i=0; i<threads; i++) {
		thrd_join(jobs[i], NULL);
	}

	/* Free resources */
	seqfclose(file);
	free(jobs);
	free(jobarg);

	return counter;
}


static KatssCounter *
count_file(const char *filename, unsigned int kmer, const char filetype)
{
	KatssCounter *counter = NULL;

	/* Open file and prepare counter & hasher */
	SeqFile read_file = seqfopen(filename, "b");
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
		still_reading = seqfread_unlocked(read_file, buffer, BUFFER_SIZE);
		buffer[still_reading] = '\0';

		katss_set_seq(hasher, buffer, filetype);
		while(katss_get_fh(hasher, &hash_value, filetype)) {
			katss_increment(counter, hash_value);
		}
	} while(still_reading == BUFFER_SIZE);

	/* If error was encountered while reading report and return NULL */
	if(still_reading == 0 && seqferrno) {
		katss_free_counter(counter);
		counter = NULL;
		error_message("katss: %d: %s", seqferrno, seqfstrerror(seqferrno));
	}

cleanup_hasher:
	free(hasher);
cleanup_file:
	seqfclose(read_file);
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
	while(seqfread(args->seqfile, buffer, BUFFER_SIZE)) {
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

/*==============================================================================
 Ushuffle counting functions
==============================================================================*/
KatssCounter *
katss_count_kmers_ushuffle(const char *filename, unsigned int kmer, int klet)
{
	/* Determine file type, or return NULL on error */
	char filetype = determine_filetype(filename);
	if(filetype == 'e' || filetype == 'N')
		return NULL;

	if(klet < 1)
		return NULL;

	KatssCounter *counter = NULL;

	/* Open file and prepare counter & hasher */
	char *buffer = s_calloc(BUFFER_SIZE, sizeof *buffer);
	char *shuf   = s_calloc(BUFFER_SIZE, sizeof *shuf);
	uint32_t hash_value;

	/* Open file */
	filetype = filetype == 'r' ? 's' : filetype;
	buffer[0] = filetype;
	SeqFile read_file = seqfopen(filename, buffer);
	if(read_file == NULL)
		goto exit;

	KatssHasher *hasher = katss_init_hasher(kmer, filetype);
	if(hasher == NULL)
		goto cleanup_file;

	counter = katss_init_counter(kmer);
	if(counter == NULL)
		goto cleanup_hasher;

	srand(1); // reset rand seed for shuffle
	while(seqfgets_unlocked(read_file, buffer, BUFFER_SIZE)) {
		int seqlen = strlen(buffer);
		shuffle(buffer, shuf, seqlen, klet);
		shuf[seqlen] = '\0';
		katss_set_seq(hasher, shuf, 'r');
		while(katss_get_fh(hasher, &hash_value, 'r')) {
			katss_increment(counter, hash_value);
		}
	}

	/* If error was encountered while reading report and return NULL */
	if(seqferrno) {
		katss_free_counter(counter);
		counter = NULL;
		error_message("katss: %d: %s", seqferrno, seqfstrerror(seqferrno));
	}

cleanup_hasher:
	free(hasher);
cleanup_file:
	seqfclose(read_file);
exit:
	free(buffer);
	free(shuf);
	return counter;
}

KatssCounter *
katss_count_kmers_ushuffle_bootstrap(const char *filename, unsigned int kmer,
                                     int klet, int sample, unsigned int *seed)
{
	/* sample should be between 1-100000 */
	sample = MAX2(sample, 1);
	sample = MIN2(sample, 100000);

	/* If not subsampling, just do regular ushuffle */
	if(sample == 100000)
		return katss_count_kmers_ushuffle(filename, kmer, klet);

	/* Check klet */
	if(klet < 1)
		return NULL;

	char filetype = determine_filetype(filename);
	if(filetype == 'e' || filetype == 'N')
		return NULL;

	KatssCounter *counter = NULL;

	/* Initialize buffer */
	char *buffer = s_calloc(BUFFER_SIZE, sizeof *buffer);
	char *shuf   = s_calloc(BUFFER_SIZE, sizeof *shuf);
	uint32_t hash_value;

	/* Open file and hasher */
	buffer[0] = filetype == 'r' ? 's' : filetype;
	SeqFile read_file = seqfopen(filename, buffer);
	if(read_file == NULL)
		goto exit;
	
	KatssHasher *hasher = katss_init_hasher(kmer, filetype);
	if(hasher == NULL)
		goto cleanup_file;
	
	counter = katss_init_counter(kmer);
	if(counter == NULL)
		goto cleanup_hasher;

	/* int to subsample from rand() */
	unsigned int local_seed;
	if(seed == NULL) {
		local_seed = time(NULL);
		seed = &local_seed;
	}

	srand(1); // reset rand seed for shuffle
	while(seqfgets_unlocked(read_file, buffer, BUFFER_SIZE)) {
		/* Pick random sequences */
		if(rand_r(seed) % 100000 >= sample)
			continue;
		/* Shuffle sequences */
		int seqlen = strlen(buffer);
		shuffle(buffer, shuf, strlen(buffer), klet);
		shuf[seqlen] = '\0'; // add null terminator since shuffle uses strncpy
		katss_set_seq(hasher, shuf, 'r');
		while(katss_get_fh(hasher, &hash_value, 'r')) {
			katss_increment(counter, hash_value);
		}
	}

	if(seqferrno) {
		error_message("katss: sample: %s\n", seqfstrerror_r(seqferrno, buffer, BUFFER_SIZE));
		katss_free_counter(counter);
		counter = NULL;
	}

cleanup_hasher:
	free(hasher);
cleanup_file:
	seqfclose(read_file);
exit:
	free(buffer);
	free(shuf);
	return counter;
}

/*==============================================================
|  Helper Functions                                            |
==============================================================*/
static char
determine_filetype(const char *file)
{
	/* Open the SeqFile, return 'e' upon error */
	SeqFile reads_file = seqfopen(file, "b");
	if(reads_file == NULL) {
		error_message("katss: %s: %s", file, seqfstrerror(seqferrno));
		seqfclose(reads_file);
		return 'N';
	}

	char buffer[BUFFER_SIZE];
	int lines_read = 0;
	int fastq_score_lines = 0;
	int fasta_score_lines = 0;
	int sequence_lines = 0;

	while (seqfgets(reads_file, buffer, BUFFER_SIZE) != NULL && lines_read < 10) {
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
    seqfclose(reads_file);

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
