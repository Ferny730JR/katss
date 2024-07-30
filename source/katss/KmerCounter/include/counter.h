#ifndef KATSS_COUNTER_H
#define KATSS_COUNTER_H

#include <stdint.h>
#include <stdatomic.h>

#ifdef __cplusplus
extern "C" {
#endif 

typedef struct KatssCounter KatssCounter;

typedef enum KATSS_TYPE {
	KATSS_INT8,
	KATSS_UINT8,
	KATSS_INT16,
	KATSS_UINT16,
	KATSS_INT32,
	KATSS_UINT32,
	KATSS_INT64,
	KATSS_UINT64,
	KATSS_FLOAT,
	KATSS_DOUBLE,
} KATSS_TYPE;


/**
 * @brief Initialize the hash tables for katss k-mer counting
 * 
 * @param kmer The size of k-mer value to count.
 * @return KatssCounter* 
 */
KatssCounter *
katss_init_counter(unsigned int kmer);


/**
 * @brief Increments the count of the hashed k-mer by one. Use the hash provided by `KmerHasher`
 * struct in `hash_function.h`
 * 
 * @param counter Pointer to KatssCounter struct
 * @param hash    Hash value of a k-mer sequence
 */
void
katss_increment(KatssCounter *counter, uint32_t hash);


/**
 * @brief Increments the count of all hash values in an array. Use the hash provided by `KmerHasher`
 * struct in `hash_functions.h`
 * 
 * @param counter     Pointer to KatssCounter struct
 * @param hash_values Array of hash values
 * @param num_values  Number of hash values in array
 */
void
katss_increments(KatssCounter *counter, uint32_t *hash_values, size_t num_values);

/**
 * @brief Decrements the count of the hashed k-mer by one. Use the hash provided by `KmerHasher`
 * struct in `hash_function.h`
 * 
 * @param counter Pointer to KatssCounter struct
 * @param hash    Hash value of a k-mer sequence
 */
void
katss_decrement(KatssCounter *counter, uint32_t hash);


/**
 * @brief Get the value associated with a key.
 * 
 * @param counter Pointer to KatssCounter struct
 * @param value   Pointer to store counts for `key`
 * @param key     K-mer to obtain counts from
 * @return int: `0` if successfully set value. `1` if `key` contained an unrecognizable char. `2`
 * if `key` is longer than the specified counter k-mer length.
 * 
 * 
 * @example
 * unsigned long count_aa;
 * katss_get(mycounter, &count_aa, "AA");
 * printf("%s: %llu\n", "AA", count_aa);
 */
int
katss_get(KatssCounter *counter, KATSS_TYPE numeric_type, void *value, const char *key);


/**
 * @brief Get the value associated with a hash value.
 * 
 * @param counter      Pointer to KatssCounter struct
 * @param numeric_type Type of `value`
 * @param value        Numeric pointer used to store count
 * @param hash         Hash to obtain count from
 * @return int `0` if successfully set value. `1` if `hash` is not in counter.
 */
int
katss_get_from_hash(KatssCounter *counter, KATSS_TYPE numeric_type, void *value, uint32_t hash);


/**
 * @brief Frees all allocated resources from KatssCounter.
 * 
 * @param counter Pointer to initialized KatssCounter struct
 */
void
katss_free_counter(KatssCounter *counter);


/**
 * @brief Count all forward-strand k-mers in a file. Currently supports fasta, fastq, and reads
 * files.
 * 
 * @param filename Name of the file containing the reads
 * @param kmer     Size of k-mers to count
 * @return KatssCounter pointer
 */
KatssCounter *katss_count_kmers(const char *filename, unsigned int kmer);


/**
 * @brief Count all forward-strand k-mers in a file. Currently supports fasta, fastq, and reads
 * files.
 * 
 * @param filename Name of the file containing the reads
 * @param kmer     Size of k-mers to counts
 * @param threads  Number of threads to use
 */
KatssCounter *katss_count_kmers_mt(const char *filename, unsigned int kmer, int threads);


/**
 * @brief Uncount a kmer
 * 
 * @param counter 
 * @param filename 
 * @param kmer 
 * @return int 
 */
int katss_uncount_kmer(KatssCounter *counter, const char *filename, const char *kmer);


/**
 * @brief Uncount a kmer multithreaded
 * 
 * @param counter 
 * @param filename 
 * @param kmer 
 * @param threads 
 * @return int 
 */
int katss_uncount_kmer_mt(KatssCounter *counter, const char *filename, const char *kmer, int threads);

#ifdef __cplusplus
}
#endif

#endif // KATSS_COUNTER_H
