#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

#include "bpp_tables.h"
#include "memory_utils.h"
#include "string_utils.h"

typedef struct Hash {
	unsigned int hash;
	unsigned int errno;
} Hash;


static Entry *
create_entry(unsigned int key, unsigned int col);

static Hash
hash(const char *key);

static void
free_entry(Entry *entry);


kmerHashTable *
init_kmer_table(unsigned int kmer, unsigned int cols)
{
	unsigned long table_size = 1 << (2 * kmer); // 4^kmer
	kmerHashTable *hash_table = s_malloc(sizeof *hash_table);

	hash_table->capacity = table_size;
	hash_table->cols = cols;
	hash_table->kmer = kmer;
	hash_table->entries = s_malloc(table_size * sizeof(Entry*));
	pthread_mutex_init(&hash_table->lock, NULL);

	for (unsigned long i = 0; i < table_size; i++) {
		hash_table->entries[i] = NULL;
	}

	return hash_table;
}


kmerHashTable *
init_bpp_table(unsigned int kmer)
{
	return init_kmer_table(kmer, kmer+1);
}


double *
kmer_get(kmerHashTable *hash_table, const char *key)
{
	Hash hash_value = hash(key);
	return hash_table->entries[hash_value.hash]->values;
}


void
kmer_add_value(kmerHashTable   *hash_table, 
                    const char      *key, 
                    double          value, 
                    unsigned int    value_index)
{
	if(value_index >= hash_table->cols) {
		error_message("value_index '%d' is greater than length of value array, which is '%d'.",
		value_index, hash_table->cols);
		exit(EXIT_FAILURE);
	}

	Hash hash_value = hash(key);
	if(hash_value.errno == 1) {
		return;
	}

	pthread_mutex_lock(&hash_table->lock);
	if(hash_table->entries[hash_value.hash] == NULL) { // make new entry if not initialized
		Entry *new_item = create_entry(hash_value.hash, hash_table->cols);
		hash_table->entries[hash_value.hash] = new_item;
	}

	hash_table->entries[hash_value.hash]->values[value_index] += value;
	pthread_mutex_unlock(&hash_table->lock);
}


static Entry *
create_entry(unsigned int key, unsigned int col)
{
	Entry *entry = s_malloc(sizeof *entry);
	entry->values = s_calloc(col, sizeof *entry->values);
	entry->num_values = col;
	entry->hash = key;

	return entry;
}


static Hash
hash(const char *key)
{
	// Assuming the key is composed of 'A', 'U/T', 'C', 'G' characters
	unsigned int hash_value = 0;
	while (*key) {
		switch(*key) {
			case 'A': hash_value = hash_value * 4;     break;
			case 'C': hash_value = hash_value * 4 + 1; break;
			case 'G': hash_value = hash_value * 4 + 2; break;
			case 'T': hash_value = hash_value * 4 + 3; break;
			case 'U': hash_value = hash_value * 4 + 3; break;
			default : return (Hash){.hash = 0, .errno = 1};
		}
		key++;
	}
	return (Hash){.hash = hash_value, .errno = 0};
}


static void
free_entry(Entry *entry)
{
	// free(entry->key);
	free(entry->values);
	free(entry);
}


void
free_kmer_table(kmerHashTable *hash_table)
{
	for (size_t i = 0; i < hash_table->capacity; i++) {
		if(hash_table->entries[i]) {
			free_entry(hash_table->entries[i]);
		}
	}
	pthread_mutex_destroy(&hash_table->lock);
	free(hash_table->entries);
	free(hash_table);
}


static void
unhash(char *key, unsigned int hash_value, int length, int is_t) 
{
	key[length] = '\0'; // Null-terminate the string

	for (int i = length - 1; i >= 0; i--) {
		switch (hash_value % 4) {
			case 0: key[i] = 'A'; break;
			case 1: key[i] = 'C'; break;
			case 2: key[i] = 'G'; break;
			case 3: key[i] = is_t ? 'T' : 'U'; break;
		}
		hash_value /= 4;
	}
}


static void
print_table_to_file(kmerHashTable *table, FILE *table_file, char sep)
{
	for(unsigned long i = 0; i < table->capacity; i++) {
		if(table->entries[i] == NULL) {
			continue;
		}

		char *key = s_malloc(table->kmer+1);
		unhash(key, table->entries[i]->hash, table->kmer, 0);
		fprintf(table_file, "%s", key);
		for(unsigned long j = 0; j < table->cols; j++) {
			fprintf(table_file, "%c%9.6f", sep, table->entries[i]->values[j]);
		}
		fprintf(table_file,"\n");
	}
}


void
kmerHashTable_to_file(kmerHashTable *table, char *name, char file_delimiter)
{
	char *filename = concat(name, ".dsv");

	FILE *table_file = fopen(filename, "w");
	if (table_file == NULL) {
		error_message("Could not write to file '%s'\n",filename);
	}

	print_table_to_file(table, table_file, file_delimiter);

	free(filename);
	fclose(table_file);
}


void
print_kmer_table(kmerHashTable *hash_table)
{
	printf("--- BEGIN KMER TABLE ---\n");
	for(unsigned long i = 0; i < hash_table->capacity; i++) {
		if(hash_table->entries[i] == NULL) {
			continue;
		}

		// printf("Key: %s, Values: ", hash_table->entries[i]->key);
		for(unsigned long j = 0; j < hash_table->cols; j++) {
			printf("%f ",hash_table->entries[i]->values[j]);
		}
		printf("\n");
	}
	printf("---- END KMER TABLE ----\n");
}
