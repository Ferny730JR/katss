#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>

#include "memory_utils.h"

#if KATSS_VERBOSE == 1
static void _error_message(const char *format);
static void _warning_message(const char *format);
#endif

void
error_message(const char *format, ...)
{
#if KATSS_VERBOSE == 1
	va_list args;
	char buffer[MAX_MESSAGE_LENGTH];

	va_start(args, format);
	vsnprintf(buffer, sizeof(buffer), format, args);
	va_end(args);

	_error_message(buffer);
}

static void 
_error_message(const char *message)
{
	char full_message[MAX_MESSAGE_LENGTH];

	snprintf(full_message, sizeof(full_message), 
	 ANSI_COLOR_RED "ERROR: " ANSI_COLOR_RESET ANSI_COLOR_BRIGHT "%s" ANSI_COLOR_RESET "\n", message);

	fprintf(stderr, "%s", full_message);
#else
	(void)format;
	return;
#endif // KATSS_VERBOSE == 1
}


void
warning_message(const char *format, ...)
{
#if KATSS_VERBOSE == 1
	va_list args;
	char buffer[MAX_MESSAGE_LENGTH];

	va_start(args, format);
	vsnprintf(buffer, sizeof(buffer), format, args);
	va_end(args);

	_warning_message(buffer);
}

static void 
_warning_message(const char *message)
{
    char full_message[MAX_MESSAGE_LENGTH];

    snprintf(full_message, sizeof(full_message), 
     ANSI_COLOR_YELLOW "WARNING: " ANSI_COLOR_RESET ANSI_COLOR_BRIGHT "%s" ANSI_COLOR_RESET "\n", message);

    fprintf(stderr, "%s", full_message);
#else
	(void)format;
	return;
#endif // KATSS_VERBOSE==1
}


void *s_malloc(size_t mem_size) {
	void *pointer = malloc(mem_size);

	if(!pointer && mem_size) {
		error_message("Could not allocate memory.");
		exit(EXIT_FAILURE);
	}

	return pointer;
}


void *s_calloc(size_t count, size_t mem_size) {
	void *pointer = calloc(count, mem_size);

	if(!pointer && mem_size && count) {
		error_message("Could not allocate memory.");
		exit(EXIT_FAILURE);
	}

	return pointer;
}


void *s_realloc(void *ptr, size_t mem_size) {
	void *new_ptr = realloc(ptr, mem_size);
	if(new_ptr == NULL) {
		error_message("Failed to reallocate memory");
		exit(EXIT_FAILURE);
	}

	return new_ptr;
}
