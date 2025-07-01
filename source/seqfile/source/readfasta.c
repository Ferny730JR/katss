/* readfasta.c - Functions for parsing fasta files
 * 
 * Copyright (c) Francisco F. Cavazos 2024
 * Subject to the MIT License
 * 
 * This file should not be used in applications. It is used to implement the
 * seqf library and is subject to change.
 */

#include "seqf_read.h"

size_t
seqf_aread(seqf_statep state, unsigned char *buffer, size_t bufsize)
{
	/* Fill buffer with data */
	size_t buffer_end;
	if((buffer_end = seqf_fill(state, buffer, --bufsize)) == 0)
		return 0;

	/* Trim sequence that was not fully read */
	if(buffer_end == bufsize) {
		while(--buffer_end && buffer[buffer_end] != '>');
		size_t offset = bufsize - buffer_end;
		if(offset > state->out_bufsiz) {
			seqferrno_ = 5;
			return 0;
		}
		memcpy(state->out_buf, buffer+buffer_end, offset);
		state->next = state->out_buf;
		state->have = offset;
	} else {
		state->have = 0;
		memset(state->out_buf, 0, state->out_bufsiz);
	}

	buffer[buffer_end] = 0;
	return buffer_end;
}

size_t
seqfaread(SeqFile file, char *buffer, size_t bufsize)
{
	seqf_statep state = (seqf_statep)file;

	mtx_lock(&state->mutex);
	size_t bytes_read = seqf_aread(state, (unsigned char *)buffer, bufsize);
	mtx_unlock(&state->mutex);

	return bytes_read;
}

size_t
seqfaread_unlocked(SeqFile file, char *buffer, size_t bufsize)
{
	return seqf_aread((seqf_statep)file, (unsigned char *)buffer, bufsize);
}

char *
seqfagets_unlocked(SeqFile file, char *buffer, size_t bufsize)
{
	/* Sanity checks & early return in case we're at eof */
	if(file == NULL)
		return NULL;
	seqf_statep state = (seqf_statep)file;
	if(state->eof)
		return NULL;

	/* Skip past fasta header info, we want the sequence */
	if(seqf_skipheader(state, '>') == NULL)
		return NULL;
	
	/* Variables to be used in obtaining the sequence */
	unsigned char *buf = (unsigned char *)buffer;
	unsigned char *eol;
	size_t left = bufsize - 1; // -1 to make space for null terminator

	/* Main loop: while there are still bytes we can fill within the buffer,
	   fill it with the sequence. If internal buffer runs out, fetch more bytes.
	   Repeat this until we reach the end of the sequence (determined by finding
	   the '>' character; the next fasta record), or we run out of available
	   bytes within the buffer. */
	if(left) do {
		if(state->have == 0 && seqf_fetch(state) != 0)
			return NULL;
		if(state->have == 0)
			break;
		if(*state->next == '>')
			break;

		seqf_shiftandcopy(state, buf, left, eol);
	} while(left);
	buf[0] = '\0';

	return buffer;
}

char *
seqfagets(SeqFile file, char *buffer, size_t bufsize)
{
	seqf_statep state = (seqf_statep)file;

	mtx_lock(&state->mutex);
	char *ret = seqfagets_unlocked(file, buffer, bufsize);
	mtx_unlock(&state->mutex);

	return ret;
}

int
seqfagetnt_unlocked(SeqFile file)
{
	seqf_statep state = (seqf_statep)file;
	if(state == NULL || state->eof)
		return EOF;
	if(state->have == 0 && seqf_fetch(state) != 0)
		return EOF;
	if(state->have == 0) /* Fetched no bytes, return EOF */
		return EOF;
	
	/* Skip past newline characters, we're interested in what comes next */
	if(*state->next == '\n') {
		state->have--;
		state->next++;
	}

	/* If start of fasta header, skip it to get to nt, or return on error */
	if(*state->next == '>' && seqf_skipline(state) == NULL)
		return EOF;

	state->have--;
	return *state->next++;
}

int
seqfagetnt(SeqFile file)
{
	seqf_statep state = (seqf_statep)file;

	mtx_lock(&state->mutex);
	int ret = seqfagetnt_unlocked(file);
	mtx_unlock(&state->mutex);

	return ret;
}
