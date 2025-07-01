/* seqf_read.c - Common internal read functions used by other seq readers
 * 
 * Copyright (c) 2024-2025 Francisco F. Cavazos
 * Subject to the MIT License
 */

#include "seqf_read.h"

#ifdef _WIN32
    #include <windows.h>
    #include <io.h>
	#include <basetsd.h>
    #define read _read
	typedef SSIZE_T ssize_t;
#else
    #include <unistd.h>
#endif


/**
 * @brief Load the state's buffer for PLAIN compression. Fills `buffer` with at 
 * most `bufsize` bytes. Number of bytes in buffer is specified by `nread`.
 * 
 * Read continuously until either the buffer is full, or an error/EOF is reached.
 * Multiple read calls are necessary since read is not guaranteed to fill the
 * buffer with the requested number of bytes. As such, keep track of how many
 * bytes it has read, and request the bytes that it needs.
 * 
 * EOF if set in the state when `read` function returns 0 (signifying EOF in 
 * file descriptor) **AND** `nread` is 0, meaning the buffer is empty.
 * 
 * @param state    File state to read from
 * @param buffer   Buffer to fill with bytes
 * @param bufsize  Number of bytes to read
 * @param nread    Number of bytes actually read
 * @return int 0 on success, -1 on error
 */
static int
seqf_loadp(seqf_statep state, unsigned char *buffer, size_t bufsize, size_t *nread)
{
	size_t left = bufsize;
	ssize_t n;
	*nread = 0;
	if(left) do {
		n = read(state->fd, buffer, left);
		if(n <= 0)
			break;
		left -= n;
	} while(left);
	if(n == -1) {
		seqferrno_ = 1;
		return -1;
	}
	*nread = bufsize - left;
	if(n == 0 && *nread == 0)
		state->eof = true;
	return 0;
}

extern int
seqf_load(seqf_statep state, unsigned char *buffer, size_t bufsize, size_t *nread)
{
	/* Process plain file */
	if(state->compression == PLAIN) {
		if(seqf_loadp(state, buffer, bufsize, nread) != 0)
			return -1;
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
			size_t nread;
			if(seqf_loadp(state, state->in_buf, state->in_bufsiz, &nread) != 0)
				return -1;
			if(nread == 0)
				break;
			state->stream.avail_in = nread;
			state->stream.next_in = state->in_buf;
		}

		/* Decompress input buffer into output */
#if defined _IGZIP_H
		ret = isal_inflate(&state->stream);
		if(ret != ISAL_DECOMP_OK && ret != ISAL_END_INPUT)
			return 3;
		left = state->stream.avail_out;
	} while(left && ret != ISAL_END_INPUT);
#else
		ret = inflate(&state->stream, Z_NO_FLUSH);
		if(ret != Z_BUF_ERROR && ret != Z_OK && ret != Z_STREAM_END) {
			seqferrno_ = 1;
			return 3;
		}
		left = state->stream.avail_out;
	} while(left && ret != Z_STREAM_END);
#endif
	*nread = bufsize - left;

	return 0;
}

extern int
seqf_fetch(seqf_statep state)
{
	if(seqf_load(state, state->out_buf, state->out_bufsiz, &state->have) != 0)
		return 1;
	state->next = state->out_buf;
	return 0;
}

extern size_t
seqf_fill(seqf_statep state, unsigned char *buffer, size_t bufsize)
{
	if(bufsize == (size_t)0)
		return (size_t)0;

	/* Declare variables */
	size_t left = bufsize;

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
	if(seqf_load(state, buffer, left, &got) != 0)
		return 0;
	left -= got;

	return bufsize - left;
}

extern unsigned char *
seqf_skipheader(seqf_statep state, char skip)
{
	size_t n;
	unsigned char *end;
	char find = skip;
	char found_skp = false;
	char found_eol = false;
	do {
		if(state->have == 0 && seqf_fetch(state) != 0)
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

extern unsigned char *
seqf_skipline(seqf_statep state)
{
	size_t n;
	unsigned char *eol;
	do {
		/* Check if bytes available in internal buffer, fetch if none */
		if(state->have == 0 && seqf_fetch(state) != 0)
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
