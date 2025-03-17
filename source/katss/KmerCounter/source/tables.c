#include <stdlib.h>
#include <float.h>
#include <stdatomic.h>
#include <limits.h>

#include "katss_core.h"
#include "counter.h"
#include "memory_utils.h"

/* Function declarations */
static void init_small_table(KatssCounter *counter, unsigned int kmer);
static void init_medium_table(KatssCounter *counter, unsigned int kmer);

/*===================================
|  Main functions (used in header)  |
===================================*/
KatssCounter *
katss_init_counter(unsigned int kmer)
{
	KatssCounter *counter = s_malloc(sizeof *counter);
	counter->kmer = kmer;
	counter->total = 0;
	// atomic_init(&counter->total, 0);
	mtx_init(&counter->lock, mtx_plain);
	counter->removed = NULL;

	if(kmer == 0 || kmer > 16) {
		error_message("KatssCounter currently does not support kmer value of '%d'.\n"
		              "Currently supported: 1-16.", kmer);
		free(counter);
		return NULL;
	}

	if(kmer <= 12)
		init_small_table(counter, kmer);
	else if(kmer <= 16)
		init_medium_table(counter, kmer);

	return counter;
}


void
katss_free_counter(KatssCounter *counter)
{
	if(counter == NULL)
		return;

	if(counter->kmer <= 12) {
		free(counter->table.small);
	} else if(counter->kmer <= 16) {
		free(counter->table.medium);
	} else {
		error_message("Kmer value of '%d' is greater than allowed range.", counter->kmer);
	}

	mtx_destroy(&counter->lock);

	katss_str_node_t *head = counter->removed;
	while(head != NULL) {
		katss_str_node_t *tmp = head;
		head = head->next;
		if(tmp->str) free(tmp->str);
		free(tmp);
	}

	free(counter);
}


void
katss_increments(KatssCounter *counter, uint32_t *hash_values, size_t num_values)
{
	mtx_lock(&counter->lock);

	if(counter->kmer <= 12)
		for(size_t i=0; i<num_values; i++)
			counter->table.small[hash_values[i]]++;
	else
		for(size_t i=0; i<num_values; i++)
			counter->table.medium[hash_values[i]]++;
	counter->total++;

	mtx_unlock(&counter->lock);
}


void
katss_increment(KatssCounter *counter, uint32_t hash)
{
	if(counter->kmer <= 12) {
		// atomic_fetch_add_explicit(&counter->table.small[hash], 1, memory_order_relaxed);
		counter->table.small[hash]++;
	} else {
		// atomic_fetch_add_explicit(&counter->table.medium[hash], 1, memory_order_relaxed);
		counter->table.medium[hash]++;
	}
	counter->total++;
	// atomic_fetch_add_explicit(&counter->total, 1, memory_order_relaxed);
}


void
katss_decrement(KatssCounter *counter, uint32_t hash)
{
	mtx_lock(&counter->lock);
	if(counter->kmer <=12) {
		counter->table.small[hash]--;
		// atomic_fetch_sub_explicit(&counter->table.small[hash], 1, memory_order_relaxed);
	} else {
		counter->table.medium[hash]--;
		// atomic_fetch_sub_explicit(&counter->table.medium[hash], 1, memory_order_relaxed);
	}
	counter->total--;
	// atomic_fetch_sub_explicit(&counter->total, 1, memory_order_relaxed);

	mtx_unlock(&counter->lock);
}


int
katss_get(KatssCounter *counter, KATSS_TYPE numeric_type, void *value, const char *key)
{
	uint32_t hash = 0, keylen = 0;

	/* Get the hash of key */
	while(*key) {
		keylen++;
		switch(*key++) {
		case 'A': hash = hash * 4;     break;
		case 'C': hash = hash * 4 + 1; break;
		case 'G': hash = hash * 4 + 2; break;
		case 'T': hash = hash * 4 + 3; break;
		case 'U': hash = hash * 4 + 3; break;
		default: return 1; /* can't hash key */
		}
	}

	/* Not a valid key if not same len as kmer */
	if(keylen != counter->kmer) {
		return 2;
	}

	/* Get value from key */

	uint64_t count;
	switch(numeric_type) {
	case KATSS_INT8:
		count = (counter->kmer <= 12 ? counter->table.small[hash] : counter->table.medium[hash]);
		*((int8_t *)value) = (int8_t)(count > INT8_MAX) ? INT8_MAX : count;
		break;
	case KATSS_UINT8:
		count = (counter->kmer <= 12 ? counter->table.small[hash] : counter->table.medium[hash]);
		*((uint8_t *)value) = (uint8_t)(count > UINT8_MAX) ? UINT8_MAX : count;
		break;
	case KATSS_INT16:
		count = (counter->kmer <= 12 ? counter->table.small[hash] : counter->table.medium[hash]);
		*((int16_t *)value) = (int16_t)(count > INT16_MAX) ? INT16_MAX : count;
		break;
	case KATSS_UINT16:
		count = (counter->kmer <= 12 ? counter->table.small[hash] : counter->table.medium[hash]);
		*((uint16_t *)value) = (uint16_t)(count > UINT16_MAX) ? UINT16_MAX : count;
		break;
	case KATSS_INT32:
		count = (counter->kmer <= 12 ? counter->table.small[hash] : counter->table.medium[hash]);
		*((int32_t *)value) = (int32_t)(count > INT32_MAX) ? INT32_MAX : count;
		break;
	case KATSS_UINT32:
		count = (counter->kmer <= 12 ? counter->table.small[hash] : counter->table.medium[hash]);
		*((uint32_t *)value) = (uint32_t)(count > UINT32_MAX) ? UINT32_MAX : count;
		break;
	case KATSS_INT64:
		count = (counter->kmer <= 12 ? counter->table.small[hash] : counter->table.medium[hash]);
		*((int64_t *)value) = (int64_t)(count > INT64_MAX) ? INT64_MAX : count;
		break;
	case KATSS_UINT64:
		count = (counter->kmer <= 12 ? counter->table.small[hash] : counter->table.medium[hash]);
		*((uint64_t *)value) = count;
		break;
	case KATSS_FLOAT:
		*((float *)value) = (float)(counter->kmer <= 12 ? counter->table.small[hash] : counter->table.medium[hash]);
		break;
	case KATSS_DOUBLE:
		*((double *)value) = (double)(counter->kmer <= 12 ? counter->table.small[hash] : counter->table.medium[hash]);
		break;
	}

	return 0;
}


int
katss_get_from_hash(KatssCounter *counter, KATSS_TYPE numeric_type, void *value, uint32_t hash)
{
	/* hash is not contained within counter */
	if(hash > counter->capacity) {
		return 1;
	}

	uint64_t count;
	switch(numeric_type) {
	case KATSS_INT8:
		count = (counter->kmer <= 12 ? counter->table.small[hash] : counter->table.medium[hash]);
		*((int8_t *)value) = (int8_t)(count > INT8_MAX) ? INT8_MAX : count;
		break;
	case KATSS_UINT8:
		count = (counter->kmer <= 12 ? counter->table.small[hash] : counter->table.medium[hash]);
		*((uint8_t *)value) = (uint8_t)(count > UINT8_MAX) ? UINT8_MAX : count;
		break;
	case KATSS_INT16:
		count = (counter->kmer <= 12 ? counter->table.small[hash] : counter->table.medium[hash]);
		*((int16_t *)value) = (int16_t)(count > INT16_MAX) ? INT16_MAX : count;
		break;
	case KATSS_UINT16:
		count = (counter->kmer <= 12 ? counter->table.small[hash] : counter->table.medium[hash]);
		*((uint16_t *)value) = (uint16_t)(count > UINT16_MAX) ? UINT16_MAX : count;
		break;
	case KATSS_INT32:
		count = (counter->kmer <= 12 ? counter->table.small[hash] : counter->table.medium[hash]);
		*((int32_t *)value) = (int32_t)(count > INT32_MAX) ? INT32_MAX : count;
		break;
	case KATSS_UINT32:
		count = (counter->kmer <= 12 ? counter->table.small[hash] : counter->table.medium[hash]);
		*((uint32_t *)value) = (uint32_t)(count > UINT32_MAX) ? UINT32_MAX : count;
		break;
	case KATSS_INT64:
		count = (counter->kmer <= 12 ? counter->table.small[hash] : counter->table.medium[hash]);
		*((int64_t *)value) = (int64_t)(count > INT64_MAX) ? INT64_MAX : count;
		break;
	case KATSS_UINT64:
		count = (counter->kmer <= 12 ? counter->table.small[hash] : counter->table.medium[hash]);
		*((uint64_t *)value) = count;
		break;
	case KATSS_FLOAT:
		*((float *)value) = (float)(counter->kmer <= 12 ? counter->table.small[hash] : counter->table.medium[hash]);
		break;
	case KATSS_DOUBLE:
		*((double *)value) = (double)(counter->kmer <= 12 ? counter->table.small[hash] : counter->table.medium[hash]);
		break;
	}

	return 0;
}


/*===================================
|  Helper functions                 |
===================================*/
static void
init_small_table(KatssCounter *counter, unsigned int kmer)
{
	counter->capacity = (1ULL << 2*kmer) - 1; /* Capacity is 4^kmer */
	counter->table.small = s_calloc(((size_t)counter->capacity+1), sizeof *counter->table.small);
	// for(uint32_t i=0; i <= counter->capacity; i++) {
	// 	atomic_init(&counter->table.small[i], 0);
	// }
}


static void
init_medium_table(KatssCounter *counter, unsigned int kmer)
{
	counter->capacity = (1ULL << 2*kmer ) - 1;
	counter->table.medium = s_calloc(((size_t)counter->capacity+1), sizeof *counter->table.medium);
	// for(uint32_t i=0; i <= counter->capacity; i++) {
	// 	atomic_init(&counter->table.medium[i], 0);
	// }
}
