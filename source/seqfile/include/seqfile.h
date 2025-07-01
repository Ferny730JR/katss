/* seqfile.h - API for the SeqFile parsing library
 *
 * Copyright (c) 2024-2025 Francisco F. Cavazos
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in 
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef SEQFILES_H
#define SEQFILES_H

#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define SEQF_VERSION "0.1.0"
#define SEQF_VERSION_NUMBER 0x00001000U
#define SEQF_VERSION_MAJOR 0
#define SEQF_VERSION_MINOR 1
#define SEQF_VERSION_PATCH 0

#define SEQFBUFSIZ 8192
#ifndef EOF
#define EOF (-1)
#endif

/**
 * @brief Get the error number encountered by SeqFile
 * 
 * @return int error number
 */
int *seqfgeterrno(void);


/**
 * @brief Error number encountered by SeqFile
 */
#define seqferrno (*seqfgeterrno())


/**
 * @brief Opaque pointer for SeqFile struct
 */
typedef struct SeqFile *SeqFile;


/**
 * @brief Open a fasta/fastq/sequence file for reading. Mode is used to specify
 * which type of file is used for reading. "a" for fasta, "q" for fastq, "s"
 * for sequence file, and "b" for binary.
 * 
 * @param path Path to the file you want to open for reading
 * @param mode Type of file being opened
 * @return SeqFile 
 */
SeqFile seqfopen(const char *path, const char *mode);


/**
 * @brief Open a fasta/fastq/sequence file for reading. Mode is used to specify
 * which type of file is used for reading. "a" for fasta, "q" for fastq, "s"
 * for sequence file, and "b" for binary. 
 * 
 * This takes an existing file descriptor (obtained through open/creat/dup/etc)
 * and uses it to process SeqFile. Since SeqFile is currently only supported
 * for file reading, the file descriptor must be opened with for reading. Once
 * `seqfclose` is called, the the file descriptor will be closed alongside the
 * SeqFile handle.
 * 
 * @param fd    File descriptor of the file you want to read
 * @param mode  Type of file being opened
 * @return SeqFile 
 */
SeqFile seqfdopen(int fd, const char *mode);


/**
 * @brief Close an SeqFile handle.
 * 
 * Removes all allocated resources associated with SeqFile.
 * 
 * @param file SeqFile to close
 */
int seqfclose(SeqFile file);


/**
 * @brief Rewind a SeqFile to read from the beginning of the file.
 * 
 * @param file SeqFile to rewind
 * @return int return code: 0 if success, failure otherwise
 */
int seqfrewind(SeqFile file);


/**
 * @brief Test if SeqFile is at the end of file.
 * 
 * @param file SeqFile pointer to test
 * @return true if at end of file
 * @return false if not at end of file
 */
bool seqfeof(SeqFile file);


/**
 * @brief Set the input buffer of the SeqFile handle
 * 
 * The input buffer is only ever used if the file is in a compressed file
 * format. If the file is not compressed, this will not affect the performance
 * of the seqf* functions.
 * 
 * @param file    SeqFile handle to set the buffer to
 * @param bufsize New size of the input buffer
 * @return int 0 on success, -1 when not enough memory was available
 */
int
seqfsetibuf(SeqFile file, size_t bufsize);


/**
 * @brief Set the output buffer of the SeqFile handle to be of size `bufsize`.
 * 
 * @param file    SeqFile handle to set the buffer to
 * @param bufsize New size of the output buffer
 * @return int 
 */
int
seqfsetobuf(SeqFile file, size_t bufsize);


/**
 * @brief Set both the input buffer of the SeqFile handle to be of size
 * `bufsize`, and the output buffer to be of size `2*bufsize`.
 * 
 * By default, the internal buffers used by SeqFile are of size `SEQFBUFSIZ`
 * for the input buffer, and `2*SEQFBUFSIZ` for the output buffer. This
 * function allows to use a custom buffer size instead of the default of
 * `SEQFBUFSIZ`.
 * 
 * @param file    SeqFile handle to set the buffer to
 * @param bufsize New size of the input/output buffers
 * @return int 0 on success. -1 when failed to set input buffer, -2 when failed
 * to set output buffer.
 */
int
seqfsetbuf(SeqFile file, size_t bufsize);


/**
 * @brief Return an allocated string detailing the error encountered from SeqFile
 * 
 * @param _seqferrno The seqferrno variable
 * @return char* String with error description
 * 
 * @note It is necessary to free the string returned from this function as memory is
 * allocated to store the string.
 */
const char *seqfstrerror(int _seqferrno);


/**
 * @brief Fill buffer with a string detailing the error code from _seqferrno
 * 
 * @param _seqferrno The seqferrno variable
 * @param buffer     Buffer to fill with error description
 * @param bufsize    Size of the buffer
 * @return int 0 on success, 1 when buffer was too small for the error message
 */
int seqfstrerror_r(int _seqferrno, char *buffer, size_t bufsize);


/**
 * @brief Read SeqFile sequences into buffer. Will continue reading until it can no
 * longer fit a full sequence within bufsize characters.
 * 
 * @param file    SeqFile pointer to read from
 * @param buffer  Buffer to write sequences to
 * @param bufsize Size of the buffer being passed
 * @return size_t Number of bytes read into buffer. 0 if end of file, or error encountered.
 */
size_t seqfread(SeqFile file, char *buffer, size_t bufsize);


/**
 * @brief Read SeqFile sequences into a buffer. Will continue reading until it can no
 * longer fit a full sequence within butsize characters.
 * 
 * @param file    SeqFile pointer to read from
 * @param buffer  Buffer to write sequences to
 * @param bufsize Size of the buffer being passed
 * @return size_t Number of bytes read into buffer. 0 if end of file, or error encountered.
 * 
 * @note
 * This function does not use a mutex to lock access to the SeqFile buffer. As such, it is not
 * thread safe. Only use in single-threaded applications
 */
size_t seqfread_unlocked(SeqFile file, char *buffer, size_t bufsize);


/** Undocumented read functions. See `seqfread()` for more information. */
size_t seqfqread(SeqFile file, char *buffer, size_t bufsize);
size_t seqfqread_unlocked(SeqFile file, char *buffer, size_t bufsize);
size_t seqfaread(SeqFile file, char *buffer, size_t bufsize);
size_t seqfaread_unlocked(SeqFile file, char *buffer, size_t bufsize);
size_t seqfsread(SeqFile file, char *buffer, size_t bufsize);
size_t seqfsread_unlocked(SeqFile file, char *buffer, size_t bufsize);


/**
 * @brief Read a record's sequence into `buffer`.
 * 
 * `seqfgets()` reads at most `bufsize - 1` characters and stores them into the `buffer` string. Reading
 * stops as soon as the end of the record it is currently processing has been reached or the end of file
 * has been reached. A `'\0'` is appended to `buffer` to ensure the string is null terminated.
 * 
 * @param file    SeqFile to read from
 * @param buffer  Buffer to fill
 * @param bufsize Size of the buffer being passed
 * @return char* Pointer to the record's sequence, or NULL if unable to find the next record.
 */
char *seqfgets(SeqFile file, char *buffer, size_t bufsize);


/**
 * @brief Read a record's sequence into `buffer`.
 * 
 * `seqfgets()` reads at most `bufsize - 1` characters and stores them into the `buffer` string. Reading
 * stops as soon as the end of the record it is currently processing has been reached or the end of file
 * has been reached. A `'\0'` is appended to `buffer` to ensure the string is null terminated.
 * 
 * @param file    SeqFile to read from
 * @param buffer  Buffer to fill
 * @param bufsize Size of the buffer being passed
 * @return char* Pointer to the record's sequence, or NULL if unable to find the next record.
 * 
 * @note
 * This function does not use a mutex to lock access to the SeqFile internal buffer. As such, it is not
 * thread-safe. Only use in single-threaded applications.
 */
char *seqfgets_unlocked(SeqFile file, char *buffer, size_t bufsize);


/** Undocumented gets functions. See seqfgets() for more information. */
char *seqfagets(SeqFile file, char *buffer, size_t bufsize);
char *seqfagets_unlocked(SeqFile file, char *buffer, size_t bufsize);
char *seqfsgets(SeqFile file, char *buffer, size_t bufsize);
char *seqfsgets_unlocked(SeqFile file, char *buffer, size_t bufsize);
char *seqfqgets(SeqFile file, char *buffer, size_t bufsize);
char *seqfqgets_unlocked(SeqFile file, char *buffer, size_t bufsize);


/**
 * @brief Read only one nucleotide from the SeqFile stream. 
 * 
 * `seqfgetnt` either returns the nucleotide or EOF when it encounters an error
 * or the end of file.
 * 
 * @param file SeqFile to read from
 * @return int Next nucleotide character, or EOF on error
 */
int seqfgetnt(SeqFile file);


/**
 * @brief Read only one nucleotide from the SeqFile stream. 
 * 
 * `seqfgetnt` either returns the nucleotide or EOF when it encounters an error
 * or the end of file.
 * 
 * @param file SeqFile to read from
 * @return int Next nucleotide character, or EOF on error
 * 
 * @note
 * This function does not use a mutex to lock access to the SeqFile internal 
 * buffer. As such, it is not thread-safe. Only use in single-threaded
 * applications.
 */
int seqfgetnt_unlocked(SeqFile file);


/* Undocumented getnt functions */
int seqfagetnt(SeqFile file);
int seqfagetnt_unlocked(SeqFile file);
int seqfqgetnt(SeqFile file);
int seqfqgetnt_unlocked(SeqFile file);
int seqfsgetnt(SeqFile file);
int seqfsgetnt_unlocked(SeqFile file);


/**
 * @brief Read the next available byte from the SeqFile stream
 * 
 * `seqfgetc` either returns the bytes, or EOF when it encounter an error of the
 * end of file.
 * 
 * @param file SeqFile to read from
 * @return int Next byte, or EOF on error
 */
int seqfgetc(SeqFile file);


/**
 * @brief Read the next available byte from the SeqFile stream
 * 
 * `seqfgetc` either returns the bytes, or EOF when it encounter an error of the
 * end of file.
 * 
 * @param file SeqFile to read from
 * @return int Next byte, or EOF on error
 * 
 * @note
 * This function does not use a mutex to lock access to the SeqFile internal 
 * buffer. As such, it is not thread-safe. Only use in single-threaded
 * applications.
 */
int seqfgetc_unlocked(SeqFile file);

#ifdef __cplusplus
}
#endif

#endif
