#include <stdlib.h>

#include "memory_utils.h"
#include "thread_safe_rand.h"

thread_safe_rand_t *
thread_safe_rand_init(void)
{
	thread_safe_rand_t *tsrand = s_malloc(sizeof *tsrand);
	mtx_init(&tsrand->lock, mtx_plain);
	tsrand->seed = 1; // mimic behavior when srand is not set
	return tsrand;
}

void
thread_safe_rand_free(thread_safe_rand_t *tsrand)
{
	if(tsrand == NULL)
		return;
	mtx_destroy(&tsrand->lock);
	free(tsrand);
}

void
thread_safe_srand(thread_safe_rand_t *tsrand, unsigned int seed)
{
	if(tsrand != NULL)
		tsrand->seed = seed;
}

int
thread_safe_rand(thread_safe_rand_t *tsrand)
{
	if(tsrand != NULL) {
		mtx_lock(&tsrand->lock);
		int val = rand_r(&tsrand->seed);
		mtx_unlock(&tsrand->lock);
		return val;
	}
	return -1;
}

int
thread_safe_rand_r(thread_safe_rand_t *tsrand, unsigned int *seed)
{
	if(tsrand != NULL) {
		mtx_lock(&tsrand->lock);
		int val = rand_r(seed);
		mtx_unlock(&tsrand->lock);
		return val;
	}
	return -1;
}
