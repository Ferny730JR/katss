#ifndef KATSS_CORE_H
#define KATSS_CORE_H

#include <stdint.h>
#include <stdbool.h>

#if (__STDC_NO_THREADS__)
#  include "tinycthread.h"
#else
#  include <threads.h>
#endif



/* Linked list used to store the removed k-mers */
typedef struct katss_str_node {
	char *str;
	struct katss_str_node *next;
} katss_str_node_t;


/* Internal structure for KatssCounter */
struct KatssCounter {
	unsigned int kmer;              /** Length of k-mer being counter */
	uint32_t capacity;              /** Number of k-mers in hash table (4^k) */
	uint64_t total;        /** Total k-mers counted so far */
	union {
		uint64_t *small;  /** for k<=12 use long to store large amounts */
		uint32_t *medium; /** K>12 use 32bit to save memory */
	} table;                       /** Table to store counts */
	katss_str_node_t *removed;     /** Linked list of removed kmers */
	mtx_t lock;
};


/* Internal structure for KatssHasher */
struct KatssHasher {
	unsigned char *sequence;      /** Sequence that is being processed */
	unsigned int kmer;            /** K-mer size to hash */
	bool end_of_seq;              /** If hasher has finished hashing the sequence */
	uint32_t mask;                /** 32-bit mask for specified k-mer length */
	uint32_t previous_hash;       /** Previous calculated hash. Used for rolling hash */
	bool has_previous;            /** Test if there is a previous hash */
	int endno;                    /** The state KatssHasher ended on while processing */
	int pos;                      /** The position to hash from. Used in case hashing was cut off early */
};
/*
Notes:
endno is used to signify where the hasher finished hashing.
Mainly for fastq/fasta files when buffering their file.
0 means it ended hashing a sequence.
1 means it ended hashing a non-sequence.
*/

#endif // KATSS_CORE_H
