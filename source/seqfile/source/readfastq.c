/* readfastq.c - Functions for parsing fastq files
 * 
 * Copyright (c) Francisco F. Cavazos 2024
 * Subject to the MIT License
 */

#include "seqf_read.h"

size_t
seqf_qread(seqf_statep state, unsigned char *buffer, size_t bufsize)
{
	/* Fill buffer with data, return 0 if nothing was filled */
	size_t buffer_end;
	if((buffer_end = seqf_fill(state, buffer, bufsize)) == 0)
		return 0;

	/* Trim sequence that was not fully read */
	if(buffer_end == bufsize) {
		register bool not_validated = true;
		do {
			/* If at fastq could not be validated, return error */
			if(buffer_end == 0) {
				seqferrno_ = 5;
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
		if(offset > state->out_bufsiz) {
			seqferrno_ = 5;
			return 0;
		}
		memcpy(state->out_buf, buffer+buffer_end, offset);
		state->have = offset;
		state->next = state->out_buf;
	} else {
		state->have = 0;
		memset(state->out_buf, 0, state->out_bufsiz);
	}
	buffer[buffer_end] = 0;
	return buffer_end;
}

size_t
seqfqread(SeqFile file, char *buffer, size_t bufsize)
{
	seqf_statep state = (seqf_statep)file;

	mtx_lock(&state->mutex);
	size_t bytes_read = seqf_qread(state, (unsigned char *)buffer, bufsize);
	mtx_unlock(&state->mutex);

	return bytes_read;
}

size_t
seqfqread_unlocked(SeqFile file, char *buffer, size_t bufsize)
{
	return seqf_qread((seqf_statep)file, (unsigned char *)buffer, bufsize);
}

char *
seqfqgets_unlocked(SeqFile file, char *buffer, size_t bufsize)
{
	if(file == NULL)
		return NULL;
	seqf_statep state = (seqf_statep)file;
	if(state->eof)
		return NULL;
	
	/* Find start of next sequence */
	if(seqf_skipheader(state, '@') == NULL)
		return NULL;

	/* Declare variables */
	unsigned char *buf = (unsigned char *)buffer;
	unsigned char *eol;
	size_t left = bufsize - 1;

	/* Fill buffer with fastq sequence */
	if(left) do {
		if(state->have == 0 && seqf_fetch(state) != 0)
			return NULL; // error in seqf_fetch
		if(state->have == 0)
			break;
		if(*state->next == '+')
			break;

		seqf_shiftandcopy(state, buf, left, eol);
	} while(left);

	seqf_skipline(state); /* Skip '+' line */
	seqf_skipline(state); /* Skip quality scores */

	/* Null terminate and return buffer */
	buf[0] = '\0';
	return buffer;
}

char *
seqfqgets(SeqFile file, char *buffer, size_t bufsize)
{
	seqf_statep state = (seqf_statep)file;

	mtx_lock(&state->mutex);
	char *ret = seqfqgets_unlocked(file, buffer, bufsize);
	mtx_unlock(&state->mutex);

	return ret;
}

int
seqfqgetnt_unlocked(SeqFile file)
{
	seqf_statep state = (seqf_statep)file;
	if(state == NULL || state->eof)
		return EOF;
	if(state->have == 0 && seqf_fetch(state) != 0)
		return EOF;
	if(state->have == 0) /* Fetched no bytes, return EOF */
		return EOF;
	
	/* Skip past newline character, we want to see what's next */
	if(*state->next == '\n') {
		state->have--;
		state->next++;
	}

	/* If start if fastq header, skip it to get to the nt */
	if(*state->next == '@' && seqf_skipline(state) == NULL)
		return EOF;
	
	/* If quality scores, skip it and get to the next nt */
	if(*state->next == '+' && seqf_skipheader(state, '@') == NULL)
		return EOF;
	
	/* In nt, return it and update next & num available bytes */
	state->have--;
	return *state->next++;
}

int
seqfqgetnt(SeqFile file)
{
	seqf_statep state = (seqf_statep)file;

	mtx_lock(&state->mutex);
	int ret = seqfqgetnt_unlocked(file);
	mtx_unlock(&state->mutex);

	return ret;
}
