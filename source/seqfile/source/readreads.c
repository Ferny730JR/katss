/* readreads.c - Functions for parsing *.reads files (raw sequences)
 * 
 * Copyright (c) Francisco F. Cavazos 2024
 * Subject to the MIT License
 */

#include "seqf_read.h"

size_t
seqf_sread(seqf_statep state, unsigned char *buffer, size_t bufsize)
{
	/* Fill buffer, leave space for null terminator */
	register size_t buffer_end;
	if((buffer_end = seqf_fill(state, buffer, --bufsize)) == 0)
		return 0;

	/* Trim sequence that was not fully read */
	if(buffer_end == bufsize) {
		while(--buffer_end && buffer[buffer_end] != '\n');
		size_t offset = bufsize - ++buffer_end; // increase to include \n
		if(offset > SEQF_CHUNK) {
			seqferrno_ = 5;
			return 0;
		}
		memcpy(state->out_buf, buffer+buffer_end, offset);
		state->have = offset;
		state->next = state->out_buf;
	} else {
		state->have = 0;
		memset(state->out_buf, 0, SEQF_CHUNK);
	}

	buffer[buffer_end] = '\0';
	return buffer_end;
}

size_t
seqfsread(SeqFile file, char *buffer, size_t bufsize)
{
	seqf_statep state = (seqf_statep)file;

	mtx_lock(&state->mutex);
	size_t bytes_read = seqf_sread(state, (unsigned char *)buffer, bufsize);
	mtx_unlock(&state->mutex);

	return bytes_read;
}

size_t
seqfsread_unlocked(SeqFile file, char *buffer, size_t bufsize)
{
	return seqf_sread((seqf_statep)file, (unsigned char *)buffer, bufsize);
}

char *
seqfsgets_unlocked(SeqFile file, char *buffer, size_t bufsize)
{
	if(file == NULL)
		return NULL;
	seqf_statep state = (seqf_statep)file;
	if(state->eof)
		return NULL;

	/* Declare variables */
	unsigned char *buf = (unsigned char *)buffer;
	unsigned char *eol;
	size_t left = bufsize - 1;

	/* Begin filling buffer */
	if(left) do {
		if(state->have == 0 && seqf_fetch(state) != 0)
			return NULL; // fetch encountered error
		if(state->have == 0)
			break;

		seqf_shiftandcopy(state, buf, left, eol);
	} while(left && eol == NULL);
	buf[0] = '\0';

	return buffer;
}

char *
seqfsgets(SeqFile file, char *buffer, size_t bufsize)
{
	seqf_statep state = (seqf_statep)file;

	mtx_lock(&state->mutex);
	char *ret = seqfsgets_unlocked(file, buffer, bufsize);
	mtx_unlock(&state->mutex);

	return ret;
}

int
seqfsgetnt_unlocked(SeqFile file)
{
	seqf_statep state = (seqf_statep)file;
	if(state == NULL || state->eof)
		return EOF;
	if(state->have == 0 && seqf_fetch(state) != 0)
		return EOF;
	if(state->have == 0) /* Fetched no bytes, return EOF */
		return EOF;
	
	/* Skip past newline character to get to nt */
	if(*state->next == '\n') {
		state->have--;
		state->next++;
	}

	/* Now in a nucleotide, return it */
	state->have--;
	return *state->next++;
}

int
seqfsgetnt(SeqFile file)
{
	seqf_statep state = (seqf_statep)file;

	mtx_lock(&state->mutex);
	int ret = seqfsgetnt_unlocked(file);
	mtx_unlock(&state->mutex);

	return ret;
}
