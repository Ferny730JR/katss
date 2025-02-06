#ifndef BPP_TABLES_H
#define BPP_TABLES_H

#include <stddef.h>

/**
 *  @brief Entries for kmerHashTable data structure. Stores the key value pair, where key is a 
 *  kmer and the value is a double array
*/
typedef struct {
	unsigned int num_values;
	unsigned int hash;
	double *values;
} Entry;


/**
 *  @brief Hash table data structure to store kmer's and associated information.
*/
struct kmerHashTable {
	unsigned long   capacity;
	unsigned int    cols;
	unsigned int    kmer;
	Entry           **entries;
};
typedef struct kmerHashTable kmerHashTable;


/**
 *  @brief Initializes a kmerHashTable of size kmer.
 * 
 *  kmerHashTable is used to store all possible k-mer in an RNA/DNA sequence. As such, the
 *  allocated capacity of the table will be 4^kmer to account for every possible k-mer sequence.
 *  Each k-mer will have an associated double array to store related information.
 * 
 *  Due to the nature of kmerHashTable, storing anything other than a k-mer will result in
 *  undefined behavior.
 * 
 *  @param kmer Size of k-mer that will be stored
 *  @param cols Number of elements that will be needed in the value array
 * 
 *  @return Pointer to the initialized kmerHashTable
*/
kmerHashTable *init_kmer_table(unsigned int kmer, unsigned int cols);


/**
 *  @brief Wrapper of init_kmer_table that creates a kmerHashTable of size kmer with kmer cols.
 * 
 *  @return Pointer to the initialized kmerHashTable
*/
kmerHashTable *init_bpp_table(unsigned int kmer);


/**
 *  Free all allocated memory in kmerHashTable.
 * 
 *  @param hash_table   kmerHashTable to free memory from.
*/
void free_kmer_table(kmerHashTable *hash_table);


/**
 *  @brief Get the value array of the associated key.
 * 
 *  @param hash_table   kmerHashTable to get values from.
 *  @param key          key to retrieve values from.
*/
double *kmer_get(kmerHashTable  *hash_table, 
                 const char     *key);


/**
 *  @brief Add a value to a key in kmerHashTable to the specified index.
 * 
 *  @param hash_table   kmerHashTable to be added to.
 *  @param key          Key value to add to.
 *  @param value        Value to be added.
 *  @param value_index  Index of the array value will be added to.
*/
void kmer_add_value(kmerHashTable   *hash_table, 
                    const char      *key, 
                    double          value, 
                    unsigned int    value_index);


/**
 *  @brief  Print contents of kmerHashTable to file.
 * 
 *  The file will contain a list of all entries in the kmerHashTable. The format for each entry
 *  will be key printed first, followed by all values in the double array.
 * 
 *  @param  table   kmerHashTable to read contents from.
 *  @param  name    Name of the file to write to (will become [name].dsv).
 *  @param  file_delimiter  Delimiter used to separate values in file.
*/
void kmerHashTable_to_file(kmerHashTable *table, 
                           char *name, 
                           char file_delimiter);


/**
 *  @brief Print contents of kmerHashTable to stdout.
 * 
 *  @param hash_table   kmerHashTable to print contents from
*/
void print_kmer_table(kmerHashTable *hash_table);

#endif  // BPP_HASH_TABLE_H
