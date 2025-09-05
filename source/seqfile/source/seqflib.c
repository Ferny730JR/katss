/* seqflib.c - seqf functions for opening a SeqFile stream
 * 
 * Copyright (c) 2024-2025 Francisco F. Cavazos
 * Subject to the MIT License
 */

#include <stdlib.h>

#include <fcntl.h>
#ifdef _WIN32
    #include <windows.h>
    #include <io.h>
    #define open _open
    #define close _close
    #define read _read
    #define lseek _lseek
    #define O_RDONLY _O_RDONLY
    #define O_WRONLY _O_WRONLY
    #define O_RDWR _O_RDWR
    #define O_CREAT _O_CREAT
#else
    #include <unistd.h>
#endif

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
	state->fd = -1;
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
seqfdopen(int fd, const char *mode)
{
	seqferrno_ = 0; // no error encountered. yet.

	/* Initialize SeqFile */
	seqf_statep seq_file = malloc(sizeof *seq_file);
	if(seq_file == NULL)
		EXIT_AND_SETERR(seq_file, 6);

	/* Init deafult values */
	init_seqfstatep(seq_file);

	/* Open file and check for errors */
	seq_file->fd = fd;
	if(fd < 0)
		EXIT_AND_SETERR(seq_file, 1);

	/* Initialize mutex */
	if(mtx_init(&seq_file->mutex, mtx_plain) != thrd_success)
		EXIT_AND_SETERR(seq_file, 2);
	seq_file->mutex_is_init = true;

	/* Set input and output buffers */
	seq_file->in_buf = malloc(SEQFBUFSIZ * sizeof *seq_file->in_buf);
	if(seq_file->in_buf == NULL)
		EXIT_AND_SETERR(seq_file, 6);
	seq_file->in_bufsiz = SEQFBUFSIZ;
	seq_file->out_buf = malloc(2*SEQFBUFSIZ * sizeof *seq_file->out_buf);
	if(seq_file->out_buf == NULL)
		EXIT_AND_SETERR(seq_file, 6);
	seq_file->out_bufsiz = 2*SEQFBUFSIZ;
	seq_file->next = seq_file->out_buf;

	/* Determine type of compression, if any */
	size_t nread = 0;
	do {
		size_t n = read(seq_file->fd, seq_file->in_buf, 2);
		if(n == -1) EXIT_AND_SETERR(seq_file, 3);
		if(n == 0) break; // reached EOF before reading magic bytes
		nread += n;
	} while(nread != 2);
	if(nread < 2) {
		seq_file->compression = PLAIN;
	} if(seq_file->in_buf[0] == 0x1F && seq_file->in_buf[1] == 0x8B) {
		seq_file->compression = GZIP;
	} else if (seq_file->in_buf[0] == 0x78 && (seq_file->in_buf[1] == 0x01 || 
	  seq_file->in_buf[1] == 0x5E || seq_file->in_buf[1] == 0x9C || 
	  seq_file->in_buf[1] == 0xDA)) {
        seq_file->compression = ZLIB;
    } else {
		seq_file->compression = PLAIN;
	}
	lseek(seq_file->fd, 0, SEEK_SET);

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

SeqFile
seqfopen(const char *path, const char *mode)
{
	int flags = O_RDONLY;
#ifdef _WIN32
	flags |= O_BINARY;
#endif
	int fd = open(path, flags); // currently only support reading
	if(fd == -1) {
		seqferrno_ = 1;
		return NULL;
	}
	return seqfdopen(fd, mode);
}

int
seqfclose(SeqFile file)
{
	if(file == NULL)
		return 1;
	int return_code = 0;
	seqf_statep state = (seqf_statep)file;
	if(state->fd > 2 && close(state->fd) == -1)
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
		return -1;
	seqf_statep state = (seqf_statep)file;
	if(lseek(state->fd, 0, SEEK_SET)==-1) {
		seqferrno_ = 1;
		return -1;
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
			return -1;
		if(ret != Z_OK)
			return -1;
		state->stream.next_in = state->in_buf;
		state->stream.avail_in = 0;
	}
#endif
	return 0;
}

bool
seqfeof(SeqFile file)
{
	return ((seqf_statep)file)->eof;
}

int
seqfsetibuf(SeqFile file, size_t bufsize)
{
	if(file == NULL)
		return -1;
	seqf_statep state = (seqf_statep)file;

	unsigned char *t = realloc(state->in_buf, bufsize);
	if(t == NULL) return -1;
	state->in_buf = t;
	state->in_bufsiz = bufsize;
	return 0;
}

int
seqfsetobuf(SeqFile file, size_t bufsize)
{
	if(file == NULL)
		return -1;
	seqf_statep state = (seqf_statep)file;

	unsigned char *t = realloc(state->out_buf, bufsize);
	if(t == NULL) return -1;
	state->out_buf = t;
	state->out_bufsiz = bufsize;
	return 0;
}

int
seqfsetbuf(SeqFile file, size_t bufsize)
{
	if(seqfsetibuf(file, bufsize) != 0)
		return -1;
	if(seqfsetobuf(file, bufsize << 1) != 0)
		return -2;
	return 0;
}
