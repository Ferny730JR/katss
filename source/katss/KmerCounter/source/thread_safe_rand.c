#ifdef _WIN32
#define _CRT_RAND_S
#endif
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
#ifdef _WIN32
		unsigned int result = tsrand->seed;
		rand_s(&result);
		tsrand->seed = result;
#else
		int result = rand_r(&tsrand->seed);
#endif
		mtx_unlock(&tsrand->lock);
		return result;
	}
	return -1;
}

int
thread_safe_rand_r(thread_safe_rand_t *tsrand, unsigned int *seed)
{
	if(tsrand != NULL) {
		mtx_lock(&tsrand->lock);
#ifdef _WIN32
		rand_s(seed);
		int result = (int)(*seed);
#else
		int result = rand_r(seed);
#endif
		mtx_unlock(&tsrand->lock);
		return result;
	}
	return -1;
}
