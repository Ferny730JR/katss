/* seqf_read.c - Common internal read functions used by other seq readers
 * 
 * Copyright (c) Francisco F. Cavazos 2024
 * Subject to the MIT License
 */

#include "seqf_read.h"

extern int
seqf_load(seqf_statep state, unsigned char *buffer, size_t bufsize, size_t *read)
{
	/* Process plain file */
	if(state->compression == PLAIN) {
		*read = fread(buffer, 1, bufsize, state->file);
		if(*read == 0 && ferror(state->file)) {
			seqferrno_ = 1;
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
			state->stream.avail_in = fread(state->in_buf, 1, SEQF_CHUNK, state->file);
			if(ferror(state->file)) {
				seqferrno_ = 1;
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
	*read = bufsize - left;

	return 0;
}

extern int
seqf_fetch(seqf_statep state)
{
	if(seqf_load(state, state->out_buf, SEQF_CHUNK, &state->have) != 0)
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
