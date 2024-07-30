#ifndef KATSS_HASH_FUNCTIONS_H
#define KATSS_HASH_FUNCTIONS_H

#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Opaque struct for KatssHasher */
typedef struct KatssHasher KatssHasher;


/**
 * @brief Initialize a k-mer hasher data structure. This is used to obtain all hash values
 * of all k-mers in the provided sequences. KmerHasher is smart, and should be able to obtain
 * hash k-mers from sequence reads, fasta, and fastq files.
 * 
 * @param kmer     k-mers lengths you want to look for
 * @param filetype Type of file the new sequence was obtained from . `a` for fasta files. `q` for
 * fastq files. `r` for raw sequences.
 * @return KatssHasher* 
 */
KatssHasher *katss_init_hasher(unsigned int kmer, char filetype);


/**
 * @brief Replace the previous sequence in the k-mer hasher struct with a new one. This is
 * particularly useful when an entire sequence has been hashed, so the KmerHasher can replace
 * the already hashed sequence with a new one. The new sequence will use the previous sequence
 * as context to hash from.
 * 
 * @param hasher   Pointer to the KmerHasher struct you want to 'feed'
 * @param sequence New sequence to replace the current one in the provided KmerHasher
 * @param filetype The type of file the new sequence was obtained from. `a` for fasta files. `q`
 * for fastq files. `r` for raw sequences.
 */
void katss_set_seq(KatssHasher *hasher, char *sequence, char filetype);


/**
 * @brief Get the next 32-bit forward-strand hash value contained in the sequence. If the hasher
 * has hashed all k-mers in the sequence, then it will return false. Works with fastq, fasta, and
 * raw sequences.
 * 
 * @param hasher    KmerHasher struct that contains the sequence information.
 * @param hash      Pointer that will contain the next hash
 * @param filetype  Type of file you are hashing from
 * @return true if `hash` was set successfully
 * @return false if `hash` there are no more hashes left
 */
bool katss_get_fh(KatssHasher *hasher, uint32_t *hash, char filetype);


/**
 * @brief Stores the k-mer associated with the provided hash_values in `key`.
 * 
 * @param key         String to store the unhashed string
 * @param hash_value  hash value to unhash
 * @param kmer        K-mer length that `hash_value` belongs to
 * @param use_t       Use nucleotide `T` (if true) or `U` (if false) in the key
 */
void katss_unhash(char *key, uint32_t hash_value, unsigned int kmer, bool use_t);


/**
 * @brief Determine if the end of sequence has been reached. If this returns true, that means there
 * are no more hash values to obtain from the provided sequence. If false, you can still hash the
 * current sequence.
 * 
 * @param hasher The KmerHasher struct that contains the sequence information.
 * @return `true` if the end of the sequence has been reached.
 * @return `false` if the end of the sequence has NOT been reached.
 */
bool katss_eos(KatssHasher *hasher);

#ifdef __cplusplus
}
#endif

#endif // KATSS_HASH_FUNCTIONS_H
