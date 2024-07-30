#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

#include "katss_core.h"
#include "hash_functions.h"
#include "memory_utils.h"


static void handle_endno(KatssHasher *hasher);
static inline int indxchr(const unsigned char *sequence, const char match);

/* Forward strand rolling hash functions */
static inline uint32_t fbh_r(KatssHasher *hasher);
static inline uint32_t fbh_a(KatssHasher *hasher);
static inline uint32_t fbh_q(KatssHasher *hasher);
static inline uint32_t frh(uint32_t previous_hash, uint32_t nt_value, uint32_t mask);


/*==========  Legend:  ==========*
0: 'A', 'a'                      |
1: 'C', 'c'                      |
2: 'G', 'g'                      |
3: 'T', 'U', 't', 'u'            |
4: '\0' (null terminator)        |
5: '>' (for fasta files)         |
6: '@' (for fastq files)         |
7: '+' (for fastq files)         |
8: '\n' (for multiline fasta)    |
9: Every other character         |
================================*/
static const uint32_t base[256] = {
	4, 9, 9, 9, 9, 9, 9, 9, 9, 9, 8, 9, 9, 9, 9, 9,  //0..15
	9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,  //16..31
	9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 7, 9, 9, 9, 9,  //32..47
	9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 5, 9,  //48..63
	6, 0, 9, 1, 9, 9, 9, 2, 9, 9, 9, 9, 9, 9, 9, 9,  //64..79
	9, 9, 9, 9, 3, 3, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,  //80..95
	9, 0, 9, 1, 9, 9, 9, 2, 9, 9, 9, 9, 9, 9, 9, 9,  //96..111
	9, 9, 9, 9, 3, 3, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,  //112..127
	9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,  //128..143
	9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,  //144..159
	9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,  //160..175
	9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,  //176..191
	9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,  //192..207
	9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,  //208..223
	9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,  //224..239
	9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,  //240..255
};

/*=================================================
| Main Hasher data structure                      |
=================================================*/

KatssHasher *
katss_init_hasher(unsigned int kmer, char filetype)
{
	KatssHasher *hasher = s_malloc(sizeof *hasher);
	hasher->mask = ((1ULL << 2*kmer) - 1);
	hasher->end_of_seq = false;
	hasher->kmer = kmer;
	hasher->sequence = NULL;
	hasher->endno = 0;
	hasher->has_previous = false;
	hasher->previous_hash = 0;
	hasher->pos = 0;
	(void)filetype; // silence compiler warnings.
	return hasher;
}


void
katss_set_seq(KatssHasher *hasher, char *sequence, char filetype)
{
	hasher->sequence = (unsigned char *)sequence;
	hasher->end_of_seq = false;
	handle_endno(hasher);
	hasher->endno = 0;
	(void)filetype; // silence compiler warnings. todo: fix this
}


void
katss_unhash(char *key, uint32_t hash_value, unsigned int kmer, bool use_t) 
{
	key[kmer] = '\0'; // Null-terminate the string

	for (int i = kmer - 1; i >= 0; i--) {
		switch (hash_value % 4) {
			case 0: key[i] = 'A'; break;
			case 1: key[i] = 'C'; break;
			case 2: key[i] = 'G'; break;
			case 3: key[i] = use_t ? 'T' : 'U'; break;
		}
		hash_value /= 4;
	}
}


bool
katss_eos(KatssHasher *hasher)
{
	return hasher->end_of_seq;
}


static void
handle_endno(KatssHasher *hasher)
{
	if(hasher->endno == 1) {
		while(*hasher->sequence && *hasher->sequence != '\n') {
			++hasher->sequence;
		}
		++hasher->sequence;
	} else if(hasher->endno == 2) {
		bool first = true;
		while(*hasher->sequence && (*hasher->sequence != '\n' || first)) {
			if(*hasher->sequence == '\n') first = false;
			++hasher->sequence;
		}
		++hasher->sequence;
	}
	if(*hasher->sequence == '\0') {
		hasher->endno = 1;
		hasher->end_of_seq = true;
	}
}

/*=================================================
| Forward-Strand Hashing Functions                |
=================================================*/

/* Get the forward hash of specified filetype */
bool
katss_get_fh(KatssHasher *hasher, uint32_t *hash, char filetype)
{
	/* If int pointer is NULL, can't modify it so return false */
	if(hash == NULL) {
		return false;
	}

	/* No previous hash to build off of, get new hash from observed */
	if(!hasher->has_previous) {
		switch(filetype) {
		case 'r': *hash = fbh_r(hasher); break;
		case 'a': *hash = fbh_a(hasher); break;
		case 'q': *hash = fbh_q(hasher); break;
		default: error_message("Filetype '%c' currently not supported.", filetype); break;
		}
		hasher->previous_hash = *hash;
		if(!hasher->end_of_seq)
			hasher->has_previous = true;

		return !hasher->end_of_seq;
	}

	/* Skip newline characters. Save time & necessary for multiline fasta */
	if(filetype != 'r' && *hasher->sequence == '\n') ++hasher->sequence;

	/* Begin actual hash computations */
	uint32_t x = base[*hasher->sequence];
	if(x < 4) {
		*hash = frh(hasher->previous_hash, x, hasher->mask);
		++hasher->sequence;
	} else if(x > 4) {
		switch(filetype) {
		case 'r': *hash = fbh_r(hasher); break;
		case 'a': *hash = fbh_a(hasher); break;
		case 'q': *hash = fbh_q(hasher); break;
		default: error_message("Filetype '%c' currently not supported.", filetype); break;
		}
	} else { // x == 4
		hasher->end_of_seq = true;
	}

	hasher->previous_hash = *hash;
	return !hasher->end_of_seq;
}


/* Forward-strand base hash of read file */
static inline uint32_t
fbh_r(KatssHasher *hasher)
{
	/* Variables for looping */
	uint32_t hash = hasher->pos ? hasher->previous_hash : 0;
	unsigned int kmer = hasher->kmer;

	/* Hash each nucleotide in sequence */
	for(unsigned int i = hasher->pos; i < kmer; i++) {
		switch (base[*hasher->sequence]) {
		case 0: hash = hash * 4;     hasher->pos++; break; // A
		case 1: hash = hash * 4 + 1; hasher->pos++; break; // C
		case 2: hash = hash * 4 + 2; hasher->pos++; break; // G
		case 3: hash = hash * 4 + 3; hasher->pos++; break; // T/U
		case 4: 
			hasher->end_of_seq = true;
			hasher->has_previous = false;
			return hash; // \0
		default: 
			hasher->pos = 0;
			hash = 0;
			i = -1;
			break;
		}
		++hasher->sequence;
	}
	hasher->pos = 0;
	return hash;
}

/* Forward-strand base hash of fasta file */
static inline uint32_t
fbh_a(KatssHasher *hasher)
{
	uint32_t hash = hasher->pos ? hasher->previous_hash : 0;
	int kmer = hasher->kmer;
	int shift;

	/* Hash each nucleotide in sequence */
	for(int i = hasher->pos; i < kmer; i++) {
		switch (base[*hasher->sequence]) {
		case 0: hash = hash * 4;     hasher->pos++; break; // A
		case 1: hash = hash * 4 + 1; hasher->pos++; break; // C
		case 2: hash = hash * 4 + 2; hasher->pos++; break; // G
		case 3: hash = hash * 4 + 3; hasher->pos++; break; // T/U
		case 4: 
			hasher->end_of_seq = true;
			hasher->has_previous = false;
			return hash; // \0
		case 5: // encountered '>' (seq info, we dont want that)
			hasher->pos = 0;
			shift = indxchr(hasher->sequence, '\n');
			if(shift == -1) {
				hasher->end_of_seq = true; 
				hasher->endno = 1; 
				hasher->has_previous = false; 
				return 0; 
			}
			hasher->sequence+=shift;
			hash = 0; i = -1;
			break;
		case 8: 
			i--;
			break; // newline, just skip!
		default: 
			hash = 0;
			i = -1;
			hasher->pos = 0;
			break;
		}
		++hasher->sequence;
	}
	hasher->pos = 0;
	return hash;
}


/* Forward-strand base hash of a fastq file 
TODO: FIX SKIPPING! CURRENTLY DOES NOT PROPERLY WORKS. */
static inline uint32_t
fbh_q(KatssHasher *hasher)
{
	uint32_t hash = hasher->pos ? hasher->previous_hash : 0;
	int kmer = hasher->kmer;
	int shift;

	/* Hash each nucleotide & handle non-sequence */
	for(int i = hasher->pos; i < kmer; i++) {
		switch (base[*hasher->sequence]) {
		case 0: hash = hash * 4;     hasher->pos++; break; // A
		case 1: hash = hash * 4 + 1; hasher->pos++; break; // C
		case 2: hash = hash * 4 + 2; hasher->pos++; break; // G
		case 3: hash = hash * 4 + 3; hasher->pos++; break; // T/U
		case 4:
			hasher->end_of_seq = true;
			hasher->has_previous = false;
			return hash; // \0
		case 6: // encountered '@' (seq info, we dont want that)
			hasher->pos = 0;
			shift = indxchr(hasher->sequence, '\n');
			if(shift == -1) {
				hasher->end_of_seq = true; 
				hasher->endno = 1; 
				hasher->has_previous = false; 
				return 0; 
			}
			hasher->sequence+=shift;
			hash = 0; i = -1;
			break;
		case 7: // encountered '+', skip two lines
			hasher->pos = 0;
			shift = indxchr(hasher->sequence, '\n');
			if(shift == -1) {
				hasher->end_of_seq = true;
				hasher->endno = 2;
				hasher->has_previous = false;
				return 0;
			}
			hasher->sequence+=shift+1;
			shift = indxchr(hasher->sequence, '\n');
			if(shift == -1) {
				hasher->end_of_seq = true;
				hasher->endno = 1;
				hasher->has_previous = false;
				return 0;
			}
			hasher->sequence+=shift;
			hash = 0; i = -1;
			break;
		case 8: 
			i--;
			break; // newline, just skip!
		default:
			hash = 0;
			i = -1;
			hasher->pos = 0;
			break;
		}
		++hasher->sequence;
	}
	hasher->pos = 0;
	return hash;
}


/* Forward-strand rolling hash */
static inline uint32_t
frh(uint32_t previous_hash, uint32_t nt_value, uint32_t mask)
{
	return ((previous_hash << 2) | nt_value) & mask;
}

/*========== Helper functions ==========*/

/**
 * @brief Find the index of the first occurrence of a char in a string
 * 
 * @param sequence   string to search in
 * @param match      char to search for
 * @return `int`     Index of match, -1 if not found
 */
static inline int
indxchr(const unsigned char *sequence, const char match)
{
	for(int indx=0; *sequence; indx++) {
		if(*sequence == match) {
			return indx;
		}
		++sequence;
	}
	return -1;
}
