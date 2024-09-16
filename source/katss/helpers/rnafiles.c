#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <errno.h>
#include <stdio.h>

#if _HAS_ISA_L_ == 1
#  include <igzip_lib.h>
#else
#  include <zlib.h>
#endif

#if (__STDC_NO_THREADS__)
#  include "tinycthread.h"
#else
#  include <threads.h>
#endif

#include "rnafiles.h"

#define RNAF_CHUNK 16384
#define MIN2(A, B)      ((A) < (B) ? (A) : (B))

static _Thread_local int rnaferrno_;

typedef enum COMPRESSION_TYPE {
	GZIP,
	ZLIB,
	PLAIN
} COMPRESSION_TYPE;

struct rnaf_state {
	FILE *file;                    /** File pointer */
	COMPRESSION_TYPE compression;  /** Type of compression, if any */
	unsigned char type;            /** Type of file, e.g FASTA, FASTQ, or reads */
#if defined _IGZIP_H               /** Use isa-l if available, otherwise zlib */
	struct inflate_state stream;   /** ISA-L Decompressor */
#else
	z_stream stream;               /** ZLIB Decompressor */
	bool stream_is_init;           /** Check is stream is initialized */
#endif

	unsigned char *in_buf;         /** Input buffer*/
	unsigned char *out_buf;        /** Output buffer */
	unsigned char *next;           /** Next available byte in output buffer */
	size_t have;                   /** Numberof bytes available in next */

	mtx_t mutex;                   /** Mutex for thread safe functions */
	bool mutex_is_init;            /** Check if mutex is initialized (for rnafclose) */

	bool eof;                      /** Flag to test if at end of rnafile */
};

typedef struct rnaf_state *rnaf_statep;

static void
init_rnafstatep(rnaf_statep state)
{
	state->file = NULL;
	state->compression = PLAIN;
	state->type = 'b';
#ifndef _IGZIP_H
	state->stream_is_init = false;
#endif
	state->in_buf = NULL;
	state->out_buf = NULL;
	state->next = NULL;
	state->have = 0;
	state->mutex_is_init = false;
	state->eof = false;
}

static COMPRESSION_TYPE
determine_compression(FILE *fp)
{
    unsigned char buffer[2];

    /* Read the first two bytes of a buffer */
    if (fread(buffer, 1, 2, fp) != 2) {
        return PLAIN; // Could not read two bytes, assume PLAIN
    }

    /* Rewind the file pointer to the beginning */
    fseek(fp, 0, SEEK_SET);

    /* Check for GZIP signature */
    if (buffer[0] == 0x1F && buffer[1] == 0x8B) {
        return GZIP;
    }

    /* Check for ZLIB signature */
    if (buffer[0] == 0x78 && (buffer[1] == 0x01 || buffer[1] == 0x5E ||
	    buffer[1] == 0x9C || buffer[1] == 0xDA)) {
        return ZLIB;
    }

    /* If not GZIP or ZLIB, we assume PLAIN file */
    return PLAIN;
}

static bool
extract_mode(rnaf_statep state, const char *mode)
{
	if(mode == NULL)
		return true;

	bool type_set = false;
	do {
		switch(*mode++) {
		case 'a': 
			if(type_set) return false;
			state->type = 'a'; break; /* fasta file */
		case 'q': 
			if(type_set) return false;
			state->type = 'q'; break; /* fastq file */
		case 's':
			if(type_set) return false;
			state->type = 's'; break; /* sequences file */
		case 'b':
			if(type_set) return false;
			state->type = 'b'; break; /* binary file*/
		case '\0': return true;
		default: return false;
		}
	} while(true);
}

static void *
error_rnafopen(rnaf_statep state, int _rnaferrno)
{
	rnaferrno_ = _rnaferrno;
	rnafclose((RnaFile)state);
	return NULL;
}

RnaFile
rnafopen(const char *filename, const char *mode)
{
	rnaferrno_ = 0; // no error encountered. yet.

	/* Initialize RnaFile */
	rnaf_statep rna_file = malloc(sizeof *rna_file);
	if(rna_file == NULL)
		return error_rnafopen(rna_file, 6);

	/* Init deafult values */
	init_rnafstatep(rna_file);

	/* Open file and check for errors */
	rna_file->file = fopen(filename, "r");
	if(rna_file->file == NULL)
		return error_rnafopen(rna_file, 1);

	/* Initialize mutex */
	if(mtx_init(&rna_file->mutex, mtx_plain) != thrd_success) 
		return error_rnafopen(rna_file, 2);
	rna_file->mutex_is_init = true;

	/* Determine type of compression, if any */
	rna_file->compression = determine_compression(rna_file->file);

	/* Set input and output buffers */
	rna_file->in_buf = malloc(RNAF_CHUNK * sizeof *rna_file->in_buf);
	if(rna_file->in_buf == NULL)
		return error_rnafopen(rna_file, 6);
	rna_file->out_buf = malloc(RNAF_CHUNK * sizeof *rna_file->out_buf);
	if(rna_file->out_buf == NULL)
		return error_rnafopen(rna_file, 6);
	rna_file->next = rna_file->out_buf;

	/* Initialize decompressor */
	if(rna_file->compression != PLAIN) {
#if defined _IGZIP_H
		isal_inflate_init(&rna_file->stream);
		rna_file->stream.crc_flag = rna_file->compression == GZIP ? ISAL_GZIP : ISAL_ZLIB;
		rna_file->stream.next_in = rna_file->in_buf;
#else
		/* allocate inflate state */
		int ret = Z_ERRNO;
		rna_file->stream.zalloc   = Z_NULL;
		rna_file->stream.zfree    = Z_NULL;
		rna_file->stream.opaque   = Z_NULL;
		rna_file->stream.avail_in = 0;
		rna_file->stream.next_in  = Z_NULL;
		if(rna_file->compression == GZIP)
			ret = inflateInit2(&rna_file->stream, 16 + MAX_WBITS);
		else if(rna_file->compression == ZLIB)
			ret = inflateInit(&rna_file->stream);
		else {
			rnaferrno_ = 1;
			rnafclose((RnaFile)rna_file);
		}
		if(ret != Z_OK)
			return error_rnafopen(rna_file, 1);
		rna_file->stream.next_in = rna_file->in_buf;
		rna_file->stream_is_init = true;
#endif
	}

	if(!extract_mode(rna_file, mode))
		return error_rnafopen(rna_file, 3);

	return (RnaFile)rna_file;
}

void
rnafclose(RnaFile file)
{
	if(file == NULL)
		return;
	rnaf_statep state = (rnaf_statep)file;
	if(state->file != NULL && fclose(state->file) == EOF)
		rnaferrno_ = 1;
	if(state->mutex_is_init)
		mtx_destroy(&state->mutex);
	if(state->in_buf)
		free(state->in_buf);
	if(state->out_buf)
		free(state->out_buf);
#ifndef _IGZIP_H
	if(state->stream_is_init)
		inflateEnd(&state->stream);
#endif
	free(state);
}

bool
rnafeof(RnaFile file)
{
	return ((rnaf_statep)file)->eof;
}

int
rnafgeterrno(void)
{
	return rnaferrno_;
}

char *
rnafstrerror_r(int _rnaferrno, char *buffer, size_t bufsize)
{
	switch(_rnaferrno) {
	case 0: strncpy(buffer, "No error was encountered.", bufsize); break;
	case 1: strerror_r(errno, buffer, bufsize); break;
	case 2: strncpy(buffer,  "Mutex failed to initialize.", bufsize); break;
	case 3: strncpy(buffer,  "Invalid mode passsed.", bufsize); break;
	case 4: strncpy(buffer,  "Read failed, could not determine type of file.", bufsize); break;
	case 5: strncpy(buffer,  "Read failed, sequence is larger than input buffer.", bufsize); break;
	case 6: strncpy(buffer,  "Out of memory.", bufsize); break;
	case 7: strncpy(buffer,  "gets failed, sequence is larger than passed buffer.", bufsize); break;
	default: strncpy(buffer, "Unrecognized error message.", bufsize); break;
	}
	return buffer;
}

char *
rnafstrerror(int _rnaferrno)
{
	char *buffer = malloc(100);
	if(buffer == NULL) {
		rnaferrno_ = 6;
		return NULL;
	}

	return rnafstrerror_r(_rnaferrno, buffer, 100);
}

/*==================================================================================================
|                                         Helper Functions                                         |
==================================================================================================*/
static int
rnaf_load(rnaf_statep state, unsigned char *buffer, size_t bufsize, size_t *read)
{
	/* Process plain file */
	if(state->compression == PLAIN) {
		*read = fread(buffer, 1, bufsize, state->file);
		if(*read == 0 && ferror(state->file)) {
			rnaferrno_ = 1;
			return 1;
		}
		if(*read == 0)
			state->eof = true;
		return 0;
	}

	/* Process compressed file */
	register int ret;
	register size_t left = bufsize;
	state->stream.next_out = buffer;
	state->stream.avail_out = bufsize;

	if(left) do {
		/* Refill input buffer if empty */
		if(state->stream.avail_in == 0) {
			state->stream.avail_in = fread(state->in_buf, 1, RNAF_CHUNK, state->file);
			if(ferror(state->file)) {
				rnaferrno_ = 1;
				return 1;
			}
			if(state->stream.avail_in == 0) {
				state->eof = true;
				break;
			}
			state->stream.next_in = state->in_buf;
		}

		/* Decompress buffer input buffer into output */
#if defined _IGZIP_H
		ret = isal_inflate(&state->stream);
		if(ret != ISAL_DECOMP_OK && ret != ISAL_END_INPUT)
			return 2;
		left = state->stream.avail_out;
	} while(left && ret != ISAL_END_INPUT);
#else
		ret = inflate(&state->stream, Z_NO_FLUSH);
		if(ret != Z_BUF_ERROR && ret != Z_OK && ret != Z_STREAM_END) {
			rnaferrno_ = 1;
			return 2;
		}
		left = state->stream.avail_out;
	} while(left && ret != Z_STREAM_END);
#endif
	*read = bufsize - left;

	return 0;
}

static int
rnaf_fetch(rnaf_statep state)
{
	if(rnaf_load(state, state->out_buf, RNAF_CHUNK, &state->have) != 0)
		return 1;
	state->next = state->out_buf;
	return 0;
}

static size_t
rnaf_fill(rnaf_statep state, unsigned char *buffer, size_t bufsize)
{
	if(bufsize == (size_t)0)
		return (size_t)0;

	/* Declare variables */
	register size_t left = bufsize; // decrement bufsize to guarantee fitting null terminator

	/* Fill buffer with decompressed bytes */
	if(state->have) {
		size_t n = MIN2(left, state->have);
		memcpy(buffer, state->out_buf, n);

		/* Move pointers */
		buffer += n;
		state->next += n;
		state->have -= n;
		left -= n;
	}

	/* Fill buffer with decompressed data */
	size_t got = 0;
	if(rnaf_load(state, buffer, left, &got) != 0)
		return 0;
	left -= got;

	return bufsize - left;
}

static unsigned char *
rnaf_skipheader(rnaf_statep state, char skip)
{
	size_t n;
	unsigned char *end;
	char find = skip;
	char found_skp = false;
	char found_eol = false;
	do {
		if(state->have == 0 && rnaf_fetch(state) != 0)
			return NULL;
		if(state->have == 0)
			return NULL;
		
		/* Try and skip */
		n = state->have;
		end = memchr(state->next, find, n);
		if(end != NULL) {
			n = (size_t)(end - state->next + 1);
			if(!found_skp) {
				find = '\n';
				found_skp = true;
			} else {
				found_eol = true;
			}
		}

		state->have -= n;
		state->next += n;
	} while(!found_eol);
	return state->next;
}

static unsigned char *
rnaf_skipline(rnaf_statep state)
{
	size_t n;
	unsigned char *eol;
	do {
		/* Check if bytes available in internal buffer, fetch if none */
		if(state->have == 0 && rnaf_fetch(state) != 0)
			return NULL;
		if(state->have == 0)
			return NULL;

		/* Try and skip to the end of the line */
		n = state->have;
		eol = memchr(state->next, '\n', n);
		if(eol != NULL)
			n = (size_t)(eol - state->next + 1);

		/* Move internal pointers */
		state->have -= n;
		state->next += n;
	} while(eol == NULL);
	return state->next;
}

/*==================================================================================================
|                                        FASTQ File Parsing                                        |
==================================================================================================*/
static size_t
rnaf_qread(rnaf_statep state, unsigned char *buffer, size_t bufsize)
{
	/* Fill buffer with data, return 0 if nothing was filled */
	size_t buffer_end;
	if((buffer_end = rnaf_fill(state, buffer, bufsize)) == 0)
		return 0;

	/* Trim sequence that was not fully read */
	if(buffer_end == bufsize) {
		register bool not_validated = true;
		do {
			/* If at fastq could not be validated, return error */
			if(buffer_end == 0) {
				rnaferrno_ = 5;
				return 0;
			}

			/* Validate by searching for '@', which must be followed by a '+' */
			while(--buffer_end && buffer[buffer_end] != '@');
			register int validate = buffer_end, count = 0;
			while(validate && count < 3) if(buffer[--validate] == '\n') count++;
			if(++validate && buffer[validate] == '+') not_validated = false;
		} while(not_validated);

		/* Copy seq that wasn't fully read into internal buffer */
		size_t offset = bufsize - buffer_end;
		if(offset > RNAF_CHUNK) {
			rnaferrno_ = 5;
			return 0;
		}
		memcpy(state->out_buf, buffer+buffer_end, offset);
		state->have = offset;
		state->next = state->out_buf;
	} else {
		state->have = 0;
		memset(state->out_buf, 0, RNAF_CHUNK);
	}
	buffer[buffer_end] = 0;
	return buffer_end;
}

size_t
rnafqread(RnaFile file, char *buffer, size_t bufsize)
{
	rnaf_statep state = (rnaf_statep)file;

	mtx_lock(&state->mutex);
	size_t bytes_read = rnaf_qread(state, (unsigned char *)buffer, bufsize);
	mtx_unlock(&state->mutex);

	return bytes_read;
}

size_t
rnafqread_unlocked(RnaFile file, char *buffer, size_t bufsize)
{
	return rnaf_qread((rnaf_statep)file, (unsigned char *)buffer, bufsize);
}

static char *
rnaf_qgets(rnaf_statep state, unsigned char *buffer, size_t bufsize)
{
	if(state == NULL || buffer == NULL || bufsize == 0)
		return NULL;
	if(state->eof)
		return NULL;
	
	/* Find start of next sequence */
	if(rnaf_skipheader(state, '@') == NULL)
		return NULL;

	/* Declare variables */
	register unsigned char *eos;
	register char *str = (char *)buffer;
	register size_t n, left = bufsize - 1;

	/* Fill buffer with fastq sequence */
	if(left) do {
		if(state->have == 0 && rnaf_fetch(state) != 0)
			return NULL; // error in rnaf_fetch
		if(state->have == 0 || *state->next == '+')
			break;

		/* Look for end of line in internal buffer */
		n = MIN2(state->have, left);
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

	rnaf_skipline(state); /* Skip '+' line */
	rnaf_skipline(state); /* Skip quality scores */

	/* Null terminate and return str */
	buffer[0] = '\0';
	return str;
}

char *
rnafqgets(RnaFile file, char *buffer, size_t bufsize)
{
	rnaf_statep state = (rnaf_statep)file;

	mtx_lock(&state->mutex);
	char *ret = rnaf_qgets(state, (unsigned char *)buffer, bufsize);
	mtx_unlock(&state->mutex);

	return ret;
}

char *
rnafqgets_unlocked(RnaFile file, char *buffer, size_t bufsize)
{
	return rnaf_qgets((rnaf_statep)file, (unsigned char *)buffer, bufsize);
}

/*==================================================================================================
|                                        FASTA File Parsing                                        |
==================================================================================================*/
static size_t
rnaf_aread(rnaf_statep state, unsigned char *buffer, size_t bufsize)
{
	/* Fill buffer with data */
	size_t buffer_end;
	if((buffer_end = rnaf_fill(state, buffer, --bufsize)) == 0)
		return 0;

	/* Trim sequence that was not fully read */
	if(buffer_end == bufsize) {
		while(--buffer_end && buffer[buffer_end] != '>');
		size_t offset = bufsize - buffer_end;
		if(offset > RNAF_CHUNK) {
			rnaferrno_ = 5;
			return 0;
		}
		memcpy(state->out_buf, buffer+buffer_end, offset);
		state->next = state->out_buf;
		state->have = offset;
	} else {
		state->have = 0;
		memset(state->out_buf, 0, RNAF_CHUNK);
	}

	buffer[buffer_end] = 0;
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
		if(state->have == 0 && rnaf_fetch(state) != 0)
			return NULL; // error encountered 
		if(state->have == 0)
			return NULL; // Nothing left to read, couldn't find header

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
				if(rnaf_fetch(state) != 0)
					return NULL; // error in fetch
				if(state->have == 0)
					return NULL; // Nothing left to read, end was not found
				rd = start = state->out_buf;
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
		return NULL;

	/* Declare variables */
	register unsigned char *eos;
	register char *str = (char *)buffer;
	register size_t n, left = bufsize - 1;

	/* Fill buffer with seq */
	if(left) do {
		if(state->have == 0 && rnaf_fetch(state) != 0)
			return NULL; // error in rnaf_fetch
		if(state->have == 0 || *state->next == '>')
			break;
		
		/* Look for end of line in internal buffer */
		n = MIN2(state->have, left);
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
	/* Fill buffer, leave space for null terminator */
	register size_t buffer_end;
	if((buffer_end = rnaf_fill(state, buffer, --bufsize)) == 0)
		return 0;

	/* Trim sequence that was not fully read */
	if(buffer_end == bufsize) {
		while(--buffer_end && buffer[buffer_end] != '\n');
		size_t offset = bufsize - ++buffer_end; // increase to include \n
		if(offset > RNAF_CHUNK) {
			rnaferrno_ = 5;
			return 0;
		}
		memcpy(state->out_buf, buffer+buffer_end, offset);
		state->have = offset;
		state->next = state->out_buf;
	} else {
		state->have = 0;
		memset(state->out_buf, 0, RNAF_CHUNK);
	}

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
	if(state->eof)
		return NULL;

	/* Declare variables */
	register size_t n, left = bufsize - 1;
	register char *str = (char *)buffer;
	register unsigned char *eos;

	/* Begin filling buffer */
	if(left) do {
		if(state->have == 0 && rnaf_fetch(state) != 0)
			return NULL; // fetch encountered error
		if(state->have == 0)
			break;

		n = MIN2(state->have, left);
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
	case 'a': return rnaf_aread(state, buffer, bufsize);
	case 'q': return rnaf_qread(state, buffer, bufsize);
	case 's': return rnaf_sread(state, buffer, bufsize);
	case 'b': {
		/* If there are bytes in internal buffer, use fill */
		if(state->have)
			return rnaf_fill(state, buffer, bufsize);
		/* Otherwise, load directly into buffer */
		size_t bytes_read;
		if(rnaf_load(state, buffer, bufsize, &bytes_read) != 0)
			return 0;
		return bytes_read;
	}
	default: 
		rnaferrno_ = 4;
		return 0;
	}
}

size_t
rnafread(RnaFile file, char *buffer, size_t bufsize)
{
	rnaf_statep state = (rnaf_statep)file;

	mtx_lock(&state->mutex);
	size_t bytes_read = rnaf_read(state, (unsigned char *)buffer, bufsize);
	mtx_unlock(&state->mutex);

	return bytes_read;
}

size_t
rnafread_unlocked(RnaFile file, char *buffer, size_t bufsize)
{
	return rnaf_read((rnaf_statep)file, (unsigned char *)buffer, bufsize);
}

static char *
rnaf_line(rnaf_statep state, unsigned char *buffer, size_t bufsize)
{
	/* Sanity checks */
	if(state == NULL || buffer == NULL || bufsize == 0)
		return NULL;
	if(state->eof)
		return NULL;

	/* Declare variables */
	register size_t n, left = bufsize - 1;
	register unsigned char *str = buffer;
	register unsigned char *eol;

	/* Begin filling buffer */
	if(left) do {
		if(state->have == 0 && rnaf_fetch(state) != 0)
			return NULL; // fetch encountered error
		if(state->have == 0)
			break;

		n = MIN2(state->have, left);
		eol = (unsigned char *)memchr(state->next, '\n', n);
		if(eol != NULL)
			n = (size_t)(eol - state->next) + 1;

		/* Copy line into buffer */
		memcpy(buffer, state->next, n);
		left -= n;
		buffer += n;
		state->have -= n;
		state->next += n;
	} while(left && eol == NULL);

	buffer[0] = '\0';
	return (char *)str;
}

static char *
rnaf_gets(rnaf_statep state, unsigned char *buffer, size_t bufsize)
{
	switch(state->type) {
	case 'a': return rnaf_agets(state, buffer, bufsize);
	case 'q': return rnaf_qgets(state, buffer, bufsize);
	case 's': return rnaf_sgets(state, buffer, bufsize);
	case 'b': return rnaf_line(state, buffer, bufsize);
	default:
		rnaferrno_ = 5;
		return NULL;
	}
}

char *
rnafgets(RnaFile file, char *buffer, size_t bufsize)
{
	rnaf_statep state = (rnaf_statep)file;

	mtx_lock(&state->mutex);
	char *ret = rnaf_gets(state, (unsigned char *)buffer, bufsize);
	mtx_unlock(&state->mutex);

	return ret;
}

char *
rnafgets_unlocked(RnaFile file, char *buffer, size_t bufsize)
{
	return rnaf_gets((rnaf_statep)file, (unsigned char *)buffer, bufsize);
}
