/* seqflib.c - seqf functions for interacting with a SeqFile pointer
 * 
 * Copyright (c) Francisco F. Cavazos 2024
 * Subject to the MIT License
 */

#include <stdlib.h>

#include "seqf_core.h"

#define EXIT_AND_SETERR(state, _seqferrno) \
	do { \
		seqferrno_ = _seqferrno; \
		seqfclose((SeqFile)state); \
		return NULL; \
	} while(0)

static void
init_seqfstatep(seqf_statep state)
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

static bool
extract_mode(seqf_statep state, const char *mode)
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

SeqFile
seqfopen(const char *path, const char *mode)
{
	seqferrno_ = 0; // no error encountered. yet.

	/* Initialize SeqFile */
	seqf_statep seq_file = malloc(sizeof *seq_file);
	if(seq_file == NULL)
		EXIT_AND_SETERR(seq_file, 6);

	/* Init deafult values */
	init_seqfstatep(seq_file);

	/* Open file and check for errors */
	seq_file->file = fopen(path, "r");
	if(seq_file->file == NULL)
		EXIT_AND_SETERR(seq_file, 1);

	/* Initialize mutex */
	if(mtx_init(&seq_file->mutex, mtx_plain) != thrd_success)
		EXIT_AND_SETERR(seq_file, 2);
	seq_file->mutex_is_init = true;

	/* Set input and output buffers */
	seq_file->in_buf = malloc(SEQF_CHUNK * sizeof *seq_file->in_buf);
	if(seq_file->in_buf == NULL)
		EXIT_AND_SETERR(seq_file, 6);
	seq_file->out_buf = malloc(SEQF_CHUNK * sizeof *seq_file->out_buf);
	if(seq_file->out_buf == NULL)
		EXIT_AND_SETERR(seq_file, 6);
	seq_file->next = seq_file->out_buf;

	/* Determine type of compression, if any */
	if(fread(seq_file->in_buf, sizeof(unsigned char), 2, seq_file->file) != 2) {
		seq_file->compression = PLAIN;
	} else if(seq_file->in_buf[0] == 0x1F && seq_file->in_buf[1] == 0x8B) {
		seq_file->compression = GZIP;
	} else if (seq_file->in_buf[0] == 0x78 && (seq_file->in_buf[1] == 0x01 || 
	  seq_file->in_buf[1] == 0x5E || seq_file->in_buf[1] == 0x9C || 
	  seq_file->in_buf[1] == 0xDA)) {
        seq_file->compression = ZLIB;
    } else {
		seq_file->compression = PLAIN;
	}
	fseek(seq_file->file, 0, SEEK_SET);

	/* Initialize decompressor */
	if(seq_file->compression != PLAIN) {
#if defined _IGZIP_H
		isal_inflate_init(&seq_file->stream);
		seq_file->stream.crc_flag = seq_file->compression == GZIP ? ISAL_GZIP : ISAL_ZLIB;
		seq_file->stream.next_in = seq_file->in_buf;
#else
		/* allocate inflate state */
		int ret = Z_ERRNO;
		seq_file->stream.zalloc   = Z_NULL;
		seq_file->stream.zfree    = Z_NULL;
		seq_file->stream.opaque   = Z_NULL;
		seq_file->stream.avail_in = 0;
		seq_file->stream.next_in  = Z_NULL;
		if(seq_file->compression == GZIP)
			ret = inflateInit2(&seq_file->stream, 16 + MAX_WBITS);
		else if(seq_file->compression == ZLIB)
			ret = inflateInit(&seq_file->stream);
		else {
			seqferrno_ = 1;
			seqfclose((SeqFile)seq_file);
		}
		if(ret != Z_OK)
			EXIT_AND_SETERR(seq_file, 1);
		seq_file->stream.next_in = seq_file->in_buf;
		seq_file->stream_is_init = true;
#endif
	}

	if(!extract_mode(seq_file, mode))
		EXIT_AND_SETERR(seq_file, 3);

	return (SeqFile)seq_file;
}

int
seqfclose(SeqFile file)
{
	if(file == NULL)
		return 1;
	int return_code = 0;
	seqf_statep state = (seqf_statep)file;
	if(state->file != NULL && fclose(state->file) == EOF)
		return_code = seqferrno_ = 1;
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
	return return_code;
}

int
seqfrewind(SeqFile file)
{
	if(file == NULL)
		return 1;
	seqf_statep state = (seqf_statep)file;
	if(fseek(state->file, 0, SEEK_SET)==-1) {
		seqferrno_ = 1;
		return 1;
	}
	state->have = 0;
	state->eof = false;
#if defined _IGZIP_H
	isal_inflate_reset(&state->stream);
	state->stream.crc_flag = state->compression == GZIP ? ISAL_GZIP : ISAL_ZLIB;
	state->stream.next_in = state->in_buf;
#else
	if(state->stream_is_init) {
		int ret;
		if(state->compression == GZIP)
			ret = inflateReset2(&state->stream, 16 + MAX_WBITS);
		else if(state->compression == ZLIB)
			ret = inflateReset(&state->stream);
		else
			return 1;
		if(ret != Z_OK)
			return 1;
		state->stream.next_in = state->in_buf;
	}
#endif
	return 0;
}

bool
seqfeof(SeqFile file)
{
	return ((seqf_statep)file)->eof;
}
