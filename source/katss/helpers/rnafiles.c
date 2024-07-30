#include <string.h>
#include <stdbool.h>
#include <errno.h>
#include <zlib.h>
#include "igzip_lib.h"

#include "rnafiles.h"
#include "tinycthread.h"
#include "memory_utils.h"

static _Thread_local int rnaferrno_;

struct rnaf_state {
	gzFile file;
	unsigned char *carry;
	unsigned char *next;
	size_t carry_size;
	size_t offset;
	size_t have;
	unsigned char type;
	mtx_t mutex;
};

typedef struct rnaf_state *rnaf_statep;
static bool extract_mode(rnaf_statep state, const char *mode);

RnaFile
rnafopen(const char *filename, const char *mode)
{
	rnaferrno_ = 0; // no error encountered. yet.

	rnaf_statep rna_file = s_malloc(sizeof *rna_file);
	rna_file->file = gzopen(filename, "r");
	if(rna_file->file == NULL) {
		rnaferrno_ = 1;
		free(rna_file);
		return NULL;
	}

	rna_file->carry_size = 8192;
	rna_file->carry = s_malloc(rna_file->carry_size * sizeof(unsigned char));
	rna_file->next = rna_file->carry;
	rna_file->offset = 0;
	rna_file->have = 0;
	rna_file->type = 'N';

	int ret = mtx_init(&rna_file->mutex, mtx_plain);
	if(ret == thrd_error) {
		rnaferrno_ = 2;
		free(rna_file->carry);
		gzclose(rna_file->file);
		free(rna_file);
		return NULL;
	}

	if(!extract_mode(rna_file, mode)) {
		mtx_destroy(&rna_file->mutex);
		gzclose(rna_file->file);
		free(rna_file->carry);
		free(rna_file);
		return NULL;
	}

	return (RnaFile)rna_file;
}

void
rnafclose(RnaFile file)
{
	if(file != NULL) {
		rnaf_statep state = (rnaf_statep)file;
		mtx_destroy(&state->mutex);
		gzclose(state->file);
		free(state->carry);
		free(state);
	}
}

bool
rnafeof(RnaFile file)
{
	return (bool)gzeof(((rnaf_statep)file)->file);
}

int
rnafgeterrno(void)
{
	return rnaferrno_;
}

char *
rnafstrerror(int rnaferrno_)
{
	char *buffer = s_malloc(1000);

	switch(rnaferrno_) {
	case 0: strncpy(buffer, "No error was encountered.", 1000); break;
	case 1: strerror_r(errno, buffer, 1000); break;
	case 2: strncpy(buffer, "Mutex failed to initialize.", 1000); break;
	case 3: strncpy(buffer, "Invalid mode passsed.", 1000); break;
	case 4: strncpy(buffer, "Read failed, could not determine type of file.", 1000); break;
	default: strncpy(buffer, "Unrecognized error message.", 1000); break;
	}
	return buffer;
}

/*====================================== Extract Mode Info =======================================*/
static bool
extract_mode(rnaf_statep state, const char *mode)
{
	if(mode == NULL) {
		return true;
	}

	do {
		switch(*mode++) {
		case 'a': state->type = 'a'; /* fasta file */
		case 'q': state->type = 'q'; /* fastq file */
		case 's': state->type = 's'; /* sequences file */
		case '\0': return true;
		default: rnaferrno_ = 3; return false;
		}
	} while(true);
}

/*==================================================================================================
|                                        RNA File Functions                                        |
==================================================================================================*/
static size_t
rnaf_fetch(rnaf_statep state)
{
	state->offset = 0;
	state->have = gzread(state->file, state->carry, state->carry_size);
	state->next = state->carry;
	return state->have;
}

/*==================================================================================================
|                                        FASTQ File Parsing                                        |
==================================================================================================*/
static size_t
rnaf_qread(rnaf_statep state, unsigned char *buffer, size_t bufsize) // TODO: Replace offset with next & have
{
	register size_t offset = state->offset;
	register unsigned char *carry = state->carry;
	bufsize = bufsize - 1; // so we can fit string terminator

	memcpy(buffer, carry, offset);
	size_t bytes_read = gzread(state->file, buffer+offset, bufsize-offset);
	register size_t buffer_end = bytes_read + offset;
	if(buffer_end == bufsize) {
		register bool not_validated = true;
		do {
			while(--buffer_end && buffer[buffer_end] != '@');
			register int validate = buffer_end, count = 0;
			while(validate && count < 3) if(buffer[--validate] == '\n') count++;
			if(++validate && buffer[validate] == '+') not_validated = false;
		} while(not_validated);
		offset = bufsize - buffer_end;
		memcpy(carry, buffer+buffer_end, offset);
		carry[offset] = 0;
	} else {
		offset = 0;
		memset(carry, 0, state->carry_size);
	}
	state->offset = offset;
	buffer[buffer_end] = 0;

	return buffer_end;
}


size_t
rnafqread(RnaFile file, char *buffer, size_t bufsize)
{
	rnaf_statep state = (rnaf_statep)file;

	mtx_lock(&state->mutex);
	size_t bytesread = rnaf_qread(state, (unsigned char *)buffer, bufsize);
	mtx_unlock(&state->mutex);

	return bytesread;
}

size_t
rnafqread_unlocked(RnaFile file, char *buffer, size_t bufsize)
{
	return rnaf_qread((rnaf_statep)file, (unsigned char *)buffer, bufsize);
}

/*==================================================================================================
|                                        FASTA File Parsing                                        |
==================================================================================================*/
static size_t
rnaf_aread(rnaf_statep state, unsigned char *buffer, size_t bufsize)
{
	register size_t offset = state->offset;
	register unsigned char *carry = state->carry;
	bufsize = bufsize - 1;

	memcpy(buffer, carry, offset);
	size_t bytes_read = gzread(state->file, buffer+offset, bufsize-offset);
	register size_t buffer_end = bytes_read + offset;
	if(buffer_end == bufsize) {
		while(--buffer_end && buffer[buffer_end] != '>');
		offset = bufsize - buffer_end;
		memcpy(carry, buffer+buffer_end, offset);
		carry[offset] = '\0';
	} else {
		offset = 0;
		memset(carry, 0, state->carry_size);
	}
	state->offset = offset;
	buffer[buffer_end] = '\0';

	return buffer_end;
}

size_t
rnafaread(RnaFile file, char *buffer, size_t bufsize)
{
	rnaf_statep state = (rnaf_statep)file;

	mtx_lock(&state->mutex);
	size_t bytes_read = rnaf_aread(state, (unsigned char *)buffer, bufsize);
	mtx_unlock(&state->mutex);

	return bytes_read;
}

size_t
rnafaread_unlocked(RnaFile file, char *buffer, size_t bufsize)
{
	return rnaf_aread((rnaf_statep)file, (unsigned char *)buffer, bufsize);
}

static char *
rnaf_skipaheader(rnaf_statep state)
{
	register bool found = false, newline = false;
	register unsigned char *rd;

	do {
		if(state->have == 0 && !rnaf_fetch(state))
			return NULL;

		/* reset read */
		rd = state->next;

		/* Find the start of fasta header */
		int n = state->have;
		unsigned char *start = memchr(rd, '>', n);

		/* If header found, find end of header */
		if(start) do {
			found = true;
			n = n - (start - rd);
			unsigned char *end = memchr(start, '\n', n);
			if(end) {
				newline = true;
				end++;
				n = end - rd;
				rd = end;
			} else {
				if(!rnaf_fetch(state))
					return NULL;
				rd = start = state->carry;
				n = state->have;
			}
		} while(!newline);
		state->have -= n;
	} while(!found);

	state->next = rd;
	return (char *)rd;
}

static char *
rnaf_agets(rnaf_statep state, unsigned char *buffer, size_t bufsize)
{
	/* Skip header information, in other words find the start of next sequence */
	if(rnaf_skipaheader(state) == NULL)
		return 0;

	/* Declare variables */
	register unsigned char *eos;
	register char *str = (char *)buffer;
	register size_t n, left = bufsize - 1;

	/* Fill buffer with seq */
	if(left) do {
		if(state->have == 0 && !rnaf_fetch(state))
			break;
		if(*state->next == '>')
			break;
		
		/* Look for end of line in internal buffer */
		n = state->have > left ? left : state->have;
		eos = (unsigned char *)memchr(state->next, '\n', n);
		if(eos != NULL)
			n = (size_t)(eos - state->next);

		/* Copy to end of seq, DONT include newline */
		memcpy(buffer, state->next, n);
		left -= n;
		buffer += n;

		/* Skip past newline if found */
		if(eos != NULL)
			n++;
		state->have -= n;
		state->next += n;
	} while(left);

	buffer[0] = '\0';
	return str;
}

char *
rnafagets(RnaFile file, char *buffer, size_t bufsize)
{
	rnaf_statep state = (rnaf_statep)file;

	mtx_lock(&state->mutex);
	char *ret = rnaf_agets(state, (unsigned char *)buffer, bufsize);
	mtx_unlock(&state->mutex);

	return ret;
}

char *
rnafagets_unlocked(RnaFile file, char *buffer, size_t bufsize)
{
	return rnaf_agets((rnaf_statep)file, (unsigned char *)buffer, bufsize);
}

/*==================================================================================================
|                                        READS File Parsing                                        |
==================================================================================================*/
static size_t
rnaf_sread(rnaf_statep state, unsigned char *buffer, size_t bufsize)
{
	register size_t offset = state->offset;
	register unsigned char *carry = state->carry;
	bufsize = bufsize - 1;

	memcpy(buffer, carry, offset);
	size_t bytes_read = gzread(state->file, buffer+offset, bufsize-offset);
	register size_t buffer_end = bytes_read + offset;
	if(buffer_end == bufsize) {
		while(--buffer_end && buffer[buffer_end] != '\n');
		offset = bufsize - buffer_end;
		memcpy(carry, buffer+buffer_end, offset);
		carry[offset] = 0;
	} else {
		offset = 0;
		memset(carry, 0, state->carry_size);
	}
	state->offset = offset;
	buffer[buffer_end] = '\0';

	return buffer_end;
}

size_t
rnafsread(RnaFile file, char *buffer, size_t bufsize)
{
	rnaf_statep state = (rnaf_statep)file;

	mtx_lock(&state->mutex);
	size_t bytes_read = rnaf_sread(state, (unsigned char *)buffer, bufsize);
	mtx_unlock(&state->mutex);

	return bytes_read;
}

size_t
rnafsread_unlocked(RnaFile file, char *buffer, size_t bufsize)
{
	return rnaf_sread((rnaf_statep)file, (unsigned char *)buffer, bufsize);
}

static char *
rnaf_sgets(rnaf_statep state, unsigned char *buffer, size_t bufsize)
{
	/* Sanity checks */
	if(state == NULL || buffer == NULL || bufsize == 0)
		return NULL;
	if(gzeof(state->file))
		return NULL;

	/* Declare variables */
	register size_t n, left = bufsize - 1;
	register char *str = (char *)buffer;
	register unsigned char *eos;

	/* Begin filling buffer */
	if(left) do {
		if(state->have == 0 && !rnaf_fetch(state))
			break;
		
		n = state->have > left ? left : state->have;
		eos = (unsigned char *)memchr(state->next, '\n', n);
		if(eos != NULL)
			n = (size_t)(eos - state->next);
		
		/* Copy sequence excluding newline */
		memcpy(buffer, state->next, n);
		left -= n;
		buffer += n;

		if(eos != NULL)
			n++;
		state->have -= n;
		state->next += n;
	} while(left && eos == NULL);
	buffer[0] = '\0';

	return str;
}

char *
rnafsgets(RnaFile file, char *buffer, size_t bufsize)
{
	rnaf_statep state = (rnaf_statep)file;

	mtx_lock(&state->mutex);
	char *ret = rnaf_sgets(state, (unsigned char *)buffer, bufsize);
	mtx_unlock(&state->mutex);

	return ret;
}

char *
rnafsgets_unlocked(RnaFile file, char *buffer, size_t bufsize)
{
	return rnaf_sgets((rnaf_statep)file, (unsigned char *)buffer, bufsize);
}

/*==================================================================================================
|                                       General File Parsing                                       |
==================================================================================================*/
static size_t
rnaf_read(rnaf_statep state, unsigned char *buffer, size_t bufsize)
{
	switch(state->type) {
	case 'a': return rnaf_aread(state, buffer, bufsize); break;
	case 'q': return rnaf_qread(state, buffer, bufsize); break;
	case 's': return rnaf_sread(state, buffer, bufsize); break;
	default: 
		rnaferrno_ = 4;
		return 0;
	}
}

size_t
rnafread(RnaFile file, char *buffer, size_t bufsize)
{
	rnaf_statep state = (rnaf_statep)file;

	size_t bytes_read;
	mtx_lock(&state->mutex);
	bytes_read = rnaf_read(state, (unsigned char *)buffer, bufsize);
	mtx_unlock(&state->mutex);

	return bytes_read;
}

size_t
rnafread_unlocked(RnaFile file, char *buffer, size_t bufsize)
{
	return rnaf_read((rnaf_statep)file, (unsigned char *)buffer, bufsize);
}

static char *
rnaf_gets(rnaf_statep state, unsigned char *buffer, size_t bufsize)
{
	switch(state->type) {
	case 'a': return rnaf_agets(state, buffer, bufsize);
	case 's': return rnaf_sgets(state, buffer, bufsize);
	case 'q': // todo: implement fastq
		warning_message("rnafgets: fastq not implemented yet.");
		return 0;
	default:
		rnaferrno_ = 4;
		return 0;
	}
}

char *
rnafgets(RnaFile file, char *buffer, size_t bufsize)
{
	rnaf_statep state = (rnaf_statep)file;

	mtx_lock(&state->mutex);
	char * ret = rnaf_gets(state, (unsigned char *)buffer, bufsize);
	mtx_unlock(&state->mutex);

	return ret;
}

char *
rnafgets_unlocked(RnaFile file, char *buffer, size_t bufsize)
{
	return rnaf_gets((rnaf_statep)file, (unsigned char *)buffer, bufsize);
}
