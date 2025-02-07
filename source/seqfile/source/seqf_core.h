/* seqf_core.h - Header to access seqf internal definition and public seqf*
 * functions
 * 
 * Copyright (c) Francisco F. Cavazos 2024
 * Subject to the MIT License
 * 
 * This file should not be used in applications. It is used to implement the
 * seqf library and is subject to change.
 */

#include <stdio.h>

#if (_HAS_ISA_L_ == 1)
#  include <igzip_lib.h>
#else
#  include <zlib.h>
#endif

#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 201112L) && !defined(__STDC_NO_THREADS__)
#  include <threads.h>
#else
#  include <tinycthread.h>
#endif

#include "seqfile.h"

#define SEQF_CHUNK 8192

extern _Thread_local int seqferrno_;

typedef enum SEQF_COMPRESSION {
	GZIP,
	ZLIB,
	PLAIN
} SEQF_COMPRESSION;

struct seqf_state {
	FILE *file;                    /** File pointer */
	SEQF_COMPRESSION compression;  /** Type of compression, if any */
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

typedef struct seqf_state *seqf_statep;
