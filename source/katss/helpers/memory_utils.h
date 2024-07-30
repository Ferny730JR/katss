#ifndef MEMORY_UTILS_H
#define MEMORY_UTILS_H

#define ANSI_COLOR_RED      "\x1b[31m"
#define ANSI_COLOR_GREEN    "\x1b[32m"
#define ANSI_COLOR_YELLOW   "\x1b[33m"
#define ANSI_COLOR_BLUE     "\x1b[34m"
#define ANSI_COLOR_MAGENTA  "\x1b[35m"
#define ANSI_COLOR_CYAN     "\x1b[36m"
#define ANSI_COLOR_BRIGHT   "\x1b[1m"
#define ANSI_COLOR_RESET    "\x1b[0m"
#define MAX_MESSAGE_LENGTH 1024

#ifndef KATSS_VERBOSE
#  define KATSS_VERBOSE 1
#endif

#define BASES "ACGU"

#include <stdlib.h>

/**
 *  @brief Get the minimum of two comparable values
 */
#define MIN2(A, B)      ((A) < (B) ? (A) : (B))

/**
 *  @brief Get the maximum of two comparable values
 */
#define MAX2(A, B)      ((A) > (B) ? (A) : (B))

/**
 *  @brief Safely allocates memory
 * 
 *  @param mem_size The size of memory to be allocated in bytes
 *  @return         A pointer to the allocated memory
*/
void *s_malloc(size_t mem_size);


/**
 *  @brief Safely allocates memory and initializes it to 0
 * 
 *  @param  count       The number of elements to be initialized
 *  @param  mem_size    The size of memory to be allocated for an element in bytes
 *  @return             A pointer to the allocated memory 
*/
void *s_calloc(size_t count, size_t mem_size);


/**
 *  @brief Safely reallocates memory
 * 
 *  @param  ptr         The pointer whose memory will be reallocated
 *  @param  mem_size    The size of memory to be allocated for the pointer in bytes
 * 
 *  @return             A new pointer to the reallocated memory
*/
void *s_realloc(void *ptr, size_t mem_size);


/**
 *  @brief Print an error message to stderr
 * 
 *  @param format   The error message to be printed
 *  @param ...      Optional arguments for the format string
*/
void error_message(const char *format, ...);


/**
 *  @brief Print a warning message to stderr
 * 
 *  @param format   The warning message to be printed
 *  @param ...      Optional arguments for the format string
*/
void warning_message(const char *format, ...);

#endif
