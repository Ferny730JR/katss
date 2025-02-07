/* seqf_read.h - Header to access seqf's internal definition and read functions
 * 
 * Copyright (c) Francisco F. Cavazos 2024
 * Subject to the MIT License
 * 
 * This file should not be used in applications. It is used to implement the
 * seqf library and is subject to change.
 */

#include <string.h> // For mem* functions

#include "seqf_core.h"


/**
 * @brief Get the minimum of two comparable values
 * 
 * Used in several *read functions
 */
#define MIN2(A, B)      ((A) < (B) ? (A) : (B))


/**
 * @brief Loads `buffer' with `bufsize' decompressed bytes and store the number
 * of decompressed bytes read into `read'. Returns 0 on success, 1 on a file
 * read error, 2 when encountered end of file (no bytes left to read), and 3 
 * when there was en error decompressing the file stream.
 * 
 * If the file is not of a compressed format, it simply tries to load the buffer
 * with the requested number of bytes.
 * 
 * The function sets the seqferrno and eof flag as appropiate when encountering
 * errors.
 * 
 * @param state   internal state pointer for the SeqFile
 * @param buffer  Buffer in which decompressed bytes will be stored in
 * @param bufsize Requested number of decompressed bytes
 * @param read    Actual number of decompressed bytes read into buffer
 * @return int 
 */
extern int seqf_load(seqf_statep state, unsigned char *buffer, size_t bufsize, size_t *read);


/**
 * @brief Fills the internal output buffer with decompressed bytes. Assumes that
 * state->have is 0 since it will overwrite everything in the output buffer. If
 * state->compression is PLAIN, then it will just copy the data from the file
 * into the output buffer. Updates state->next to point to the first byte of the
 * output buffer and sets state->have to the size of the internal output buffer.
 * 
 * On success return 0, otherwise return 1.
 * 
 * @param state Internal state pointer for the SeqFile
 * @return int  
 */
extern int seqf_fetch(seqf_statep state);


/**
 * @brief Fills `buffer' with `bufsize' decompressed bytes. First, it will fill
 * buffer with the next available bytes from the internal output buffer. If not
 * enough bytes were available to fill the buffer, then it will decompress and
 * load the remaining number of requested bytes directly into the buffer.
 * 
 * @param state   Internal state pointer for the SeqFile
 * @param buffer  Buffer in which decompressed bytes will be stored in
 * @param bufsize Requested number of bytes
 * @return size_t Number of bytes read into `buffer'
 */
extern size_t seqf_fill(seqf_statep state, unsigned char *buffer, size_t bufsize);


/**
 * @brief Skip to the start of sequence data by searching for and skipping a 
 * header line.
 * 
 * This function searches within the internal buffer of a `SeqFile` state for
 * the character specified by `skip` (typically used to identify headers in
 * fasta/fastq files). It advances the reading position past this line to allow
 * direct access to the sequence data.
 * 
 * @param state Pointer to the internal `SeqFile` state
 * @param skip  Character to identify the start of a header line.
 * 
 * @return unsigned char* Pointer to the next position in the buffer after the 
 *         header line. Returns `NULL` if the header line is not found or if an
 *         error occurs.
 */
extern unsigned char *seqf_skipheader(seqf_statep state, char skip);


/**
 * @brief Skip the current line in the internal buffer of a SeqFile state.
 * 
 * This function advances the reading position within the internal buffer of the 
 * `SeqFile` state, moving past the current line to the start of the next line.
 * It is typically used to bypass unwanted lines, allowing the next read
 * operation to begin from the following line.
 * 
 * @param state Pointer to the internal `SeqFile` state
 * 
 * @return unsigned char* Pointer to the start of the next line in the buffer 
 *         after the skipped line. Returns `NULL` if the end of the buffer is 
 *         reached or an error occurs.
 */
extern unsigned char *seqf_skipline(seqf_statep state);


/**
 * @brief Function boilerplate; define a helpful macro to not have to rewrite it
 * everytime.
 * 
 * This copies the line within the state's internal buffer into buf. The amount
 * to copy is determined by `left`, the amount of bytes within the buf, or the
 * number of bytes available within the internal buffer - whichever is smaller.
 * After copying, shift the state's internal buffer and buf to point to the next
 * available empty byte.
 * 
 * @param state SeqFile internal state to modify
 * @param buf   Buffer to copy data into
 * @param left  Number of bytes left in `buf`
 * @param eol   char pointer to store end of line
 */
#define seqf_shiftandcopy(state, buf, left, eol) \
	/* Get maximum bytes we are allowed to read */ \
	size_t n = MIN2(state->have, left); \
\
	/* Find how many bytes are in the line, if eol is found */ \
	eol = memchr(state->next, '\n', n); \
	if(eol != NULL) \
		n = (size_t)(eol - state->next); \
\
	/* Copy the line (without \n) and shift internal pointer */ \
	memcpy(buf, state->next, n); \
	left -= n; \
	buf  += n; \
\
	/* Skip past newline within internal buffer if eol was found */ \
	if(eol != NULL) \
		n++; \
	state->have -= n; \
	state->next += n
// end seqfshiftcpy


/* Undocumented functions. Used for file-specific reading */

size_t seqf_qread(seqf_statep state, unsigned char *buffer, size_t bufsize);
size_t seqf_aread(seqf_statep state, unsigned char *buffer, size_t bufsize);
size_t seqf_sread(seqf_statep state, unsigned char *buffer, size_t bufsize);

char *seqf_qgets(seqf_statep state, unsigned char *buffer, size_t bufsize);
char *seqf_agets(seqf_statep state, unsigned char *buffer, size_t bufsize);
char *seqf_sgets(seqf_statep state, unsigned char *buffer, size_t bufsize);
