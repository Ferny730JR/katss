#ifndef THREAD_SAFE_RAND_H
#define THREAD_SAFE_RAND_H

#if (__STDC_NO_THREADS__)
#  include "tinycthread.h"
#else
#  include <threads.h>
#endif

struct thread_safe_rand_t {
	mtx_t lock;
	unsigned int seed;
};
typedef struct thread_safe_rand_t thread_safe_rand_t;


/**
 * @brief 
 * 
 * @return thread_safe_rand_t* 
 */
thread_safe_rand_t *thread_safe_rand_init(void);


/**
 * @brief 
 * 
 * @param tsr 
 */
void thread_safe_rand_free(thread_safe_rand_t *tsr);


/**
 * @brief Equivalent to `srand`.
 * 
 * Set the seed to be used by `thread_safe_rand()`.
 * 
 * @param tsr
 * @param seed 
 */
void thread_safe_srand(thread_safe_rand_t *tsr, unsigned int seed);


/**
 * @brief Equivalent to `rand`.
 * 
 * Return a random number
 * 
 * @param tsr thread_safe_rand_t struct
 * @return int random number
 */
int thread_safe_rand(thread_safe_rand_t *tsr);


/**
 * @brief Equivalent to `rand_r()`
 * 
 * Return a random number defined by seed
 * 
 * @param tsr   thread_safe_rand_t struct
 * @param seed  Seed to use for random number
 * @return int Random number
 */
int thread_safe_rand_r(thread_safe_rand_t *tsr, unsigned int *seed);

#endif // THREAD_SAFE_RAND_H
