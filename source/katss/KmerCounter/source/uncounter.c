#include <zlib.h>
#include <stdbool.h>
#include <errno.h>
#include <string.h>

#include "katss_core.h"
#include "counter.h"
#include "rnafiles.h"
#include "hash_functions.h"
#include "memory_utils.h"
#include "seqseq.h"

#define BUFFER_SIZE 65536

struct DecrementValues {
	int shift;
	int start;
};

struct threadinfo {
	RnaFile file;
	KatssCounter *counter;
	char *kmer;
	char *(*find)(const char *, const char *);
	size_t(*read)(RnaFile, char *, size_t);
	char *(*proc)(KatssCounter *, char *, const char *);
};

typedef struct DecrementValues DecrementValues;
typedef struct threadinfo threadinfo;

/*==================== file specific uncounting functions ====================*/
static int uncount_kmer_fasta(KatssCounter *counter, const char *filename, const char *kmer);
static int uncount_kmer_fastq(KatssCounter *counter, const char *filename, const char *kmer);
static int uncount_kmer_reads(KatssCounter *counter, const char *filename, const char *kmer);
static int remove_kmer(void *arg);

/*========================= Line processing functions =========================*/
static char *process_line_fasta(KatssCounter *counter, char *found, const char *pat);
static DecrementValues decrement_kmer_fasta(KatssCounter *counter, const char *sequence, const char *pat, int min_start, int max_end);
static char *process_line(KatssCounter *counter, char *found, const char *pat);
static DecrementValues decrement_kmer(KatssCounter *counter, const char *sequence, const char *pat, int min_start, int max_end);

/*========================== File parsing functions ==========================*/
static inline bool nhash(const char *key, uint32_t *hash_value, int start, int length);
static inline void cross_out(char *s1, const char *s2);
static inline void cross_out_fasta(char *s1, const char *s2);
static inline int subindx(const char *s1, const char *s2);
static inline int subindx_fasta(const char *s1, const char *s2);
static char determine_filetype(const char *file);
static bool is_nucleotide(char character);
static void push(KatssCounter *counter, const char *str);


/*==================================================================================================
|                                         Public Functions                                         |
==================================================================================================*/
int
katss_uncount_kmer(KatssCounter *counter, const char *filename, const char *kmer)
{
	char filetype = determine_filetype(filename);
	if(filetype == 'N') {
		return -1;
	} else if (filetype == 'e') {
		error_message("Unable to read sequence from file.\nCurrent supported file types are:"
		              " FASTA, FASTQ, and file containing sequences per line.");
		return -1;
	}

	int num_removed;
	switch (filetype) {
	case 'a':
		num_removed = uncount_kmer_fasta(counter, filename, kmer);
		break;
	case 'q':
		num_removed = uncount_kmer_fastq(counter, filename, kmer);
		break;
	case 's':
		num_removed = uncount_kmer_reads(counter, filename, kmer);
		break;
	default:
		error_message("katss_uncount_kmer: This error message should be impossible to reach. "
		              "If you received this message, please email francisco.cavazos03@gmail.com");
		return -1;
	}

	/* Add the kmer to removed */
	push(counter, kmer);
	return num_removed;
}


int
katss_uncount_kmer_mt(KatssCounter *counter, const char *filename, const char *kmer, int threads)
{
	char filetype = determine_filetype(filename);
	if(filetype == 'N') {
		return -1;
	} else if (filetype == 'e') {
		error_message("Unable to read sequence from file.\nCurrent supported file types are:"
		              " FASTA, FASTQ, and file containing sequences per line.");
		return -1;
	}

	/* Get the total before */
	int previous_total;
	katss_get(counter, KATSS_INT32, &previous_total, kmer);

	/* Open RnaFile for reading */
	char mode[2] = { 0 };
	mode[0] = filetype;
	RnaFile file = rnafopen(filename, mode);
	if(file == NULL) {
		char *error = rnafstrerror(rnaferrno);
		warning_message("rnafopen: error %d: %s",rnaferrno,error);
		free(error);
		return -1;
	}

	/* Create threads for reading */
	threads = threads < 1 ? 1 : threads;
	thrd_t *jobs = s_malloc(threads * sizeof *jobs);
	threadinfo *jobarg = s_malloc(threads * sizeof *jobarg);
	for(int i=0; i<threads; i++) {
		jobarg[i].counter = counter;
		jobarg[i].file = file;
		jobarg[i].kmer = (char *)kmer;
		switch(filetype) {
		case 'a': 
			jobarg[i].read = rnafaread;
			jobarg[i].find = seqlseqa;
			jobarg[i].proc = process_line_fasta;
			break;
		case 'q':
			jobarg[i].read = rnafqread;
			jobarg[i].find = seqlseqq;
			jobarg[i].proc = process_line;
			break;
		case 's':
			jobarg[i].read = rnafsread;
			jobarg[i].find = seqlseq;
			jobarg[i].proc = process_line;
			break;
		default:
			rnafclose(file);
			free(jobarg);
			free(jobs);
			return -1;
		}

		/* Start threads */
		thrd_create(&jobs[i], remove_kmer, &jobarg[i]);
	}

	/* Join all threads */
	for(int i=0; i<threads; i++) {
		thrd_join(jobs[i], NULL);
	}

	/* Free allocated resources */
	rnafclose(file);
	free(jobarg);
	free(jobs);

	/* Add kmer to removed list */
	push(counter, kmer);
	int current_total;
	katss_get(counter, KATSS_INT32, &current_total, kmer);
	return previous_total - current_total;
}

/*==================================================================================================
|                                        Private Functions                                         |
==================================================================================================*/
static int
uncount_kmer_fasta(KatssCounter *counter, const char *filename, const char *kmer) {
	RnaFile rnafile = rnafopen(filename, "a");
	if(rnafile == NULL) {
		char *error = rnafstrerror(rnaferrno);
		error_message("rnafopen: %s", error);
		free(error);
	}

	uint64_t previous_total = counter->total;
	char buffer[BUFFER_SIZE] = { 0 };
	while(rnafagets_unlocked(rnafile, buffer, BUFFER_SIZE)) {
		process_line(counter, buffer, kmer);
	}
	rnafclose(rnafile);

	return previous_total - counter->total;
}


static int
uncount_kmer_fastq(KatssCounter *counter, const char *filename, const char *kmer)
{
	RnaFile rnafile = rnafopen(filename, "q");
	if(rnafile == NULL) {
		char *error = rnafstrerror(rnaferrno);
		error_message("rnafopen: %s", error);
		free(error);
		return -1;
	}

	uint64_t previous_total = counter->total;
	char buffer[BUFFER_SIZE] = { 0 };
	register char *ptr;
	while(rnafqread_unlocked(rnafile, buffer, BUFFER_SIZE)) {
		ptr = buffer;
		while((ptr = seqlseqq(ptr, kmer)) != NULL) {
			ptr = process_line(counter, ptr, kmer);
		}
	}
	rnafclose(rnafile);

	return previous_total - counter->total;
}


static int
uncount_kmer_reads(KatssCounter *counter, const char *filename, const char *kmer)
{
	RnaFile rnafile = rnafopen(filename, "s");
	if(rnafile == NULL) {
		char *error = rnafstrerror(rnaferrno);
		error_message("rnafopen: %s", error);
		free(error);
		return -1;
	}

	uint64_t previous_total = counter->total;
	char buffer[BUFFER_SIZE] = { 0 };
	register char *ptr;
	while(rnafsread_unlocked(rnafile, buffer, BUFFER_SIZE)) {
		ptr = buffer;
		while((ptr = seqlseq(ptr, kmer)) != NULL) {
			ptr = process_line(counter, ptr, kmer);
		}
	}
	rnafclose(rnafile);

	return previous_total - counter->total;
}


static int
remove_kmer(void *arg)
{
	threadinfo *rec = (threadinfo *)arg;
	char *buffer = s_malloc(BUFFER_SIZE * sizeof(char));
	register char *ptr;
	while(rec->read(rec->file, buffer, BUFFER_SIZE)) {
		ptr = buffer;
		while((ptr = rec->find(ptr, rec->kmer)) != NULL) {
			ptr = rec->proc(rec->counter, ptr, rec->kmer);
		}
	}
	free(buffer);
	return 0;
}


static char *
process_line_fasta(KatssCounter *counter, char *found, const char *pat)
{
	int end = 0;
	bool found_newseq = false;
	for(end=0; found[end]; end++) {
		if(found[end] == '>') {
			found_newseq = true;
			found[end] = '\0';
			break;
		}
	}

	/* Remove sequences in line */
	katss_str_node_t *cur = counter->removed;
	while(cur != NULL) {
		cross_out_fasta(found, cur->str);
		cur = cur->next;
	}

	/* Begin uncounting pat */
	int shift = 0;
	int num_kmers_in_seq = end - counter->kmer + 1;
	DecrementValues vals = {.shift = 0, .start = 0};

	while(shift < num_kmers_in_seq) {
		vals = decrement_kmer_fasta(counter, found+shift, pat, vals.start, num_kmers_in_seq-shift);
		shift += vals.shift;
	}

	/* Undo null termination, return pointer to end of seq */
	if(found_newseq) found[end] = '>';
	return found+end;
}


static DecrementValues
decrement_kmer_fasta(KatssCounter *counter, const char *sequence, const char *pat, int min_start, int max_end)
{
	DecrementValues vals = {.shift = max_end, .start = 0};

	int pat_indx = subindx_fasta(sequence, pat);
	if(pat_indx == -1) {
		return vals;
	}

	int end = pat_indx;
	int pat_len = strlen(pat);

	for(int i=0; i<pat_len; i++) {
		if(sequence[end] == '\n') i--;
		else if(sequence[end] == '\0') break;
		end++;
	}

	end = end > max_end ? max_end : end;

	int start = pat_indx;
	for(int i=0; i < (int)counter->kmer-1; i++) {
		if(sequence[start] == '\n') i--;
		start--;
	}
	start = start < min_start ? min_start : start; 

	uint32_t hash_value;
	while(start < end) {
		if(sequence[start] == '\n') {
			start++;
			continue;
		}
		bool is_valid = nhash(sequence, &hash_value, start++, counter->kmer);
		if(!is_valid) { continue; }
		katss_decrement(counter, hash_value);
	}

	vals.shift = end;
	vals.start = end - vals.shift;

	return vals;
}

static char *
process_line(KatssCounter *counter, char *found, const char *pat)
{
	int end = 0;
	bool found_newline = false;
	for(end=0; found[end]; end++) {
		if(found[end] == '\n') {
			found_newline = true;
			found[end] = '\0';
			break;
		}
	}

	/* Remove sequences in line */
	katss_str_node_t *cur = counter->removed;
	while(cur != NULL) {
		cross_out(found, cur->str);
		cur = cur->next;
	}

	/* Begin uncounting pat */
	int shift = 0;
	int num_kmers_in_seq = end - counter->kmer + 1;
	DecrementValues vals = {.shift = 0, .start = 0};

	while(shift < num_kmers_in_seq) {
		vals = decrement_kmer(counter, found+shift, pat, vals.start, num_kmers_in_seq-shift);
		shift += vals.shift;
	}

	/* Undo null termination, return pointer to end of line */
	if(found_newline) found[end] = '\n';
	return found+end;
}


static DecrementValues
decrement_kmer(KatssCounter *counter, const char *sequence, const char *pat, int min_start, int max_end)
{
	DecrementValues vals = {.shift = max_end, .start = 0};

	int pat_indx = subindx(sequence, pat);
	if(pat_indx == -1) {
		return vals;
	}

	int pat_len = strlen(pat);
	int end = pat_indx + pat_len;
	end = end > max_end ? max_end : end;

	int start = pat_indx - counter->kmer + 1;
	start = start < min_start ? min_start : start; 

	uint32_t hash_value;
	while(start < end) {
		bool is_valid = nhash(sequence, &hash_value, start++, counter->kmer);
		if(!is_valid) { continue; }
		katss_decrement(counter, hash_value);
	}

	vals.shift = pat_indx + pat_len;
	vals.start = end - vals.shift;

	return vals;
}


static char
determine_filetype(const char *file)
{
	/* Open the gzFile, return 'N' upon error */
	gzFile reads_file = gzopen(file, "r");
	if(reads_file == NULL) {
		error_message("katss: %s: %s", file, strerror(errno));
		gzclose(reads_file);
		return 'N';
	}

	char buffer[BUFFER_SIZE];
	int lines_read = 0;
	int fastq_score_lines = 0;
	int sequence_lines = 0;

	while (gzgets(reads_file, buffer, BUFFER_SIZE) != NULL && lines_read < 10) {
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
			gzclose(reads_file);
			return 'a';
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

    gzclose(reads_file);

    if (fastq_score_lines >= 2) {
        return 'q'; // fastq file
    } else if (sequence_lines == 10) {
        return 's'; // raw sequences file
    } else {
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


static inline bool
nhash(const char *key, uint32_t *hash_value, int start, int length)
{
	register uint32_t hash_value_ = 0;
	register char *ptr = (char *)key+start;
	for(int i=0; i<length; i++) {
		switch(*ptr++) {
		case 'A': hash_value_ = hash_value_ * 4;     break;
		case 'C': hash_value_ = hash_value_ * 4 + 1; break;
		case 'G': hash_value_ = hash_value_ * 4 + 2; break;
		case 'T': hash_value_ = hash_value_ * 4 + 3; break;
		case 'U': hash_value_ = hash_value_ * 4 + 3; break;
		case '\n': i--; break;
		default : return false;
		}
	}
	*hash_value = hash_value_;
	return true;
}


static inline int
subindx(const char *s1, const char *s2) {
    char *s = seqseq(s1, s2);
    if(s) {
        return s - s1;
    } else {
        return -1;
    }
}


static inline void
cross_out(char *s1, const char *s2) {
	int s2_len = (int)strlen(s2);
    int indx;
    while((indx = subindx(s1, s2)) != -1) {
        for(int i=indx; i<indx+s2_len; i++) {
            s1[i]='X';
        }
    }
}


static inline int
subindx_fasta(const char *s1, const char *s2) {
	char *s = seqseqa(s1, s2);
	if(s) {
		return s - s1;
	} else {
		return -1;
	}
}


static inline void
cross_out_fasta(char *s1, const char *s2)
{
	int s2_len = (int)strlen(s2);
	int indx;
	while((indx = subindx_fasta(s1, s2)) != -1) {
		for(int i=indx; i<indx+s2_len; i++) {
			if(s1[i] == '\0') break;
			else if(s1[i] == '\n') indx++;
			else s1[i] = 'X';
		}
	}
}


static void
push(KatssCounter *counter, const char *str)
{
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
