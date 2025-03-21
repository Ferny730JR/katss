#include <stdbool.h>
#include <errno.h>
#include <string.h>

#if (__STDC_NO_THREADS__)
#  include "tinycthread.h"
#else
#  include <threads.h>
#endif

#include "katss_core.h"
#include "counter.h"
#include "hash_functions.h"
#include "memory_utils.h"
#include "seqfile.h"
#include "seqseq.h"
#include "ushuffle.h"

#define BUFFER_SIZE 65536U

struct threadinfo {
	SeqFile seqfile;
	KatssCounter *counter;
	char filetype;
};
typedef struct threadinfo threadinfo;

static char determine_filetype(const char *filename);
static inline void cross_out(char *s1, const char *s2, char filetype);
static void kctr_push(KatssCounter *counter, const char *str);

int
katss_recount_kmer(KatssCounter *counter, const char *filename, const char *remove)
{
	int ret = 0;
	char filetype = determine_filetype(filename);
	if(filetype == 'e' || filetype == 'N')
		return 1;

	/* Clear counter */
	uint64_t total = ((uint64_t)counter->capacity) + 1;
	if(counter->kmer <= 12)
		memset(counter->table.small,  0x00, total * sizeof(uint64_t));
	else
		memset(counter->table.medium, 0x00, total * sizeof(uint32_t));
	
	/* Push kmer to remove to counter */
	kctr_push(counter, remove);

	/* Open SeqFile for reading */
	char mode[2] = { 0 };
	mode[0] = filetype == 'r' ? 's' : filetype;
	SeqFile read_file = seqfopen(filename, mode);
	if(read_file == NULL) { /* Error opening SeqFile */
		error_message("katss: seqfopen: %s\n", seqfstrerror(seqferrno));
		return 2;
	}

	/* Initialize hasher */
	KatssHasher *hasher = katss_init_hasher(counter->kmer, filetype);
	if(hasher == NULL) {
		seqfclose(read_file);
		return 3;
	}

	char *buffer = s_malloc(BUFFER_SIZE+1);
	size_t still_reading;
	uint32_t hash_value;

	/* Begin recounting */
	do {
		still_reading = seqfread_unlocked(read_file, buffer, BUFFER_SIZE);

		/* Remove sequences in line */
		katss_str_node_t *cur = counter->removed;
		while(cur != NULL) {
			cross_out(buffer, cur->str, filetype);
			cur = cur->next;
		}

		katss_set_seq(hasher, buffer, filetype);
		while(katss_get_fh(hasher, &hash_value, filetype)) {
			katss_increment(counter, hash_value);
		}
	} while(still_reading);

	/* If error was encountered while reading report and return NULL */
	if(still_reading == 0 && seqferrno) {
		ret = 4;
		error_message("katss: %d: %s", seqferrno, seqfstrerror(seqferrno));
	}

	/* Cleanup */
	free(hasher);
	free(buffer);
	seqfclose(read_file);

	return ret;
}

int
katss_recount_kmer_shuffle(KatssCounter *counter, const char *file, int klet, const char *remove)
{
	int ret = 0;
	char filetype = determine_filetype(file);
	if(filetype == 'e' || filetype == 'N')
		return 1;

	/* Clear counter */
	uint64_t total = ((uint64_t)counter->capacity) + 1;
	if(counter->kmer <= 12)
		memset(counter->table.small,  0x00, total * sizeof(uint64_t));
	else
		memset(counter->table.medium, 0x00, total * sizeof(uint32_t));
	
	/* Push kmer to remove to counter */
	kctr_push(counter, remove);

	/* Open SeqFile for reading */
	char mode[2] = { 0 };
	mode[0] = filetype == 'r' ? 's' : filetype;
	SeqFile read_file = seqfopen(file, mode);
	if(read_file == NULL) { /* Error opening SeqFile */
		error_message("katss: seqfopen: %s\n", seqfstrerror(seqferrno));
		return 2;
	}

	/* Initialize hasher */
	KatssHasher *hasher = katss_init_hasher(counter->kmer, filetype);
	if(hasher == NULL) {
		seqfclose(read_file);
		return 3;
	}

	char *buffer = s_malloc(BUFFER_SIZE);
	char *shuf   = s_malloc(BUFFER_SIZE);
	size_t still_reading;
	uint32_t hash_value;

	/* Begin recounting */
	srand(1); // reset seed for ushuffle
	while(seqfgets_unlocked(read_file, buffer, BUFFER_SIZE)) {
		/* Shuffle the sequence */
		int seqlen = strlen(buffer);
		shuffle(buffer, shuf, seqlen, klet);
		shuf[seqlen] = '\0'; // null terminate shuf since shuffle uses strncpy

		/* Remove sequences in line */
		katss_str_node_t *cur = counter->removed;
		while(cur != NULL) {
			cross_out(buffer, cur->str, filetype);
			cur = cur->next;
		}

		/* Count the kemrs */
		katss_set_seq(hasher, shuf, filetype);
		while(katss_get_fh(hasher, &hash_value, filetype)) {
			katss_increment(counter, hash_value);
		}
	}

	/* If error was encountered while reading report and return NULL */
	if(seqferrno) {
		ret = 4;
		error_message("katss: %d: %s", seqferrno, seqfstrerror(seqferrno));
	}

	/* Cleanup */
	free(hasher);
	free(buffer);
	free(shuf);
	seqfclose(read_file);

	return ret;
}

static int
recount_mt(void *arg)
{
	threadinfo *args = (threadinfo *)arg;
	char *buffer = s_malloc(BUFFER_SIZE);

	/* Hasher to hash k-mers */
	KatssHasher *hasher = katss_init_hasher(args->counter->kmer, '\0');
	if(hasher == NULL)
		return 1;

	/* Megabyte to store counts */
	size_t num_counts = 250000;
	uint32_t *hash_values = s_malloc(num_counts * sizeof *hash_values);
	size_t cur_hash = 0;

	/* Begin re-counting */
	while(seqfread(args->seqfile, buffer, BUFFER_SIZE)) {
		/* Remove unwanted k-mers */
		katss_str_node_t *cur = args->counter->removed;
		while(cur != NULL) {
			cross_out(buffer, cur->str, args->filetype);
			cur = cur->next;
		}

		/* Count the k-mers */
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

int
katss_recount_kmer_mt(KatssCounter *counter, const char *filename, const char *remove, int threads)
{
	int ret = 0;

	/* Check type of file, or throw error if not supported */
	char filetype = determine_filetype(filename);
	if(filetype == 'e' || filetype == 'N')
		return 1;

	/* Clear counter */
	uint64_t total = ((uint64_t)counter->capacity) + 1;
	if(counter->kmer <= 12)
		memset(counter->table.small,  0x00, total * sizeof(uint64_t));
	else
		memset(counter->table.medium, 0x00, total * sizeof(uint32_t));
	
	/* Push kmer to remove to counter */
	kctr_push(counter, remove);

	/* Set minimum/maximum number of threads */
	threads = MAX2(threads, 1);
	threads = MIN2(threads, 128);

	/* Open SeqFile for reading */
	char mode[2] = { 0 };
	mode[0] = filetype == 'r' ? 's' : filetype;
	SeqFile read_file = seqfopen(filename, mode);
	if(read_file == NULL) { /* Error opening SeqFile */
		error_message("katss: seqfopen: %s\n", seqfstrerror(seqferrno));
		return 2;
	}

	/* Begin preparing threads */
	threadinfo *jobarg = s_malloc(threads * sizeof *jobarg);
	thrd_t *jobs = s_malloc(threads * sizeof *jobs);

	for(int i=0; i<threads; i++) {
		jobarg[i].seqfile = read_file;
		jobarg[i].counter = counter;
		jobarg[i].filetype = filetype;

		/* Start threads */
		thrd_create(&jobs[i], recount_mt, &jobarg[i]);
	}

	for(int i=0; i<threads; i++) {
		thrd_join(jobs[i], &ret);
	}

	/* Free resources */
	seqfclose(read_file);
	free(jobs);
	free(jobarg);

	return ret;
}

/*==================================================================================================
|                                         Helper Functions                                         |
==================================================================================================*/
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
	/* Open the SeqFile, return 'e' upon error */
	SeqFile reads_file = seqfopen(file, "b");
	if(reads_file == NULL) {
		error_message("katss: %s: %s", file, strerror(errno));
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

static inline void
cross_out(char *s1, const char *s2, char filetype) {
	register size_t s2_len = strlen(s2);
	register char *ptr = s1;
	if(filetype == 'a') {
		while((ptr = seqseqa(ptr, s2)) != NULL) {
			memset(ptr, 'X', s2_len);
		}
	} else {
		while((ptr = seqseq(ptr, s2)) != NULL) {
			memset(ptr, 'X', s2_len);
		}
	}
}

static void
kctr_push(KatssCounter *counter, const char *str)
{
	if(str == NULL)
		return;

	if(counter->removed == NULL) {
		counter->removed = s_malloc(sizeof(katss_str_node_t));
		counter->removed->next = NULL;
		counter->removed->str = strdup(str);
		return;
	}

	katss_str_node_t *cur = counter->removed;
	while(cur->next != NULL) {
		cur = cur->next;
	}

	cur->next = s_malloc(sizeof(katss_str_node_t));
	cur->next->str = strdup(str);
	cur->next->next = NULL;
}
