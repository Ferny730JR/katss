/* seqfread.c - seqf functions for general file reading of SeqFile's
 * 
 * Copyright (c) Francisco F. Cavazos 2024
 * Subject to the MIT License
 */

#include "seqf_read.h"

static size_t
seqf_read(seqf_statep state, unsigned char *buffer, size_t bufsize)
{
	switch(state->type) {
	case 'a': return seqf_aread(state, buffer, bufsize);
	case 'q': return seqf_qread(state, buffer, bufsize);
	case 's': return seqf_sread(state, buffer, bufsize);
	case 'b': return seqf_fill(state, buffer, bufsize);
	default: 
		seqferrno_ = 4;
		return 0;
	}
}

size_t
seqfread(SeqFile file, char *buffer, size_t bufsize)
{
	seqf_statep state = (seqf_statep)file;

	mtx_lock(&state->mutex);
	size_t bytes_read = seqf_read(state, (unsigned char *)buffer, bufsize);
	mtx_unlock(&state->mutex);

	return bytes_read;
}

size_t
seqfread_unlocked(SeqFile file, char *buffer, size_t bufsize)
{
	return seqf_read((seqf_statep)file, (unsigned char *)buffer, bufsize);
}

static char *
seqf_line(seqf_statep state, unsigned char *buffer, size_t bufsize)
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
		if(state->have == 0 && seqf_fetch(state) != 0)
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

char *
seqfgets_unlocked(SeqFile file, char *buffer, size_t bufsize)
{
	if(file == NULL)
		return NULL;
	seqf_statep state = (seqf_statep)file;

	switch(state->type) {
	case 'a':
		return seqfagets_unlocked(file, buffer, bufsize);
	case 'q':
		return seqfqgets_unlocked(file, buffer, bufsize);
	case 's':
		return seqfsgets_unlocked(file, buffer, bufsize);
	case 'b':
		return seqf_line(state, (unsigned char *)buffer, bufsize);
	default:
		seqferrno_ = 5;
		return NULL;
	}
}

char *
seqfgets(SeqFile file, char *buffer, size_t bufsize)
{
	if(file == NULL)
		return NULL;
	seqf_statep state = (seqf_statep)file;

	mtx_lock(&state->mutex);
	char *ret = seqfgets_unlocked(file, buffer, bufsize);
	mtx_unlock(&state->mutex);

	return ret;
}

int
seqfgetc_unlocked(SeqFile file)
{
	seqf_statep state = (seqf_statep)file;
	if(state == NULL || state->eof)
		return EOF;
	if(state->have == 0 && seqf_fetch(state) != 0)
		return EOF;
	if(state->have == 0)
		return EOF;
	state->have--;
	return *state->next++;	
}

int
seqfgetc(SeqFile file)
{
	seqf_statep state = (seqf_statep)file;

	mtx_lock(&state->mutex);
	int ret = seqfgetc_unlocked(file);
	mtx_unlock(&state->mutex);

	return ret;
}

int
seqfgetnt_unlocked(SeqFile file)
{
	switch(((seqf_statep)file)->type) {
	case 'a': return seqfagetnt_unlocked(file);
	case 'q': return seqfqgetnt_unlocked(file);
	case 's': return seqfsgetnt_unlocked(file);
	case 'b': return seqfgetc_unlocked(file);
	default:  return -1;
	}
}

int
seqfgetnt(SeqFile file)
{
	seqf_statep state = (seqf_statep)file;

	mtx_lock(&state->mutex);
	int ret = seqfgetnt_unlocked(file);
	mtx_unlock(&state->mutex);

	return ret;
}
