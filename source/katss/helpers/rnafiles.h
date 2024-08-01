#ifndef KATSS_RNAFILES_H
#define KATSS_RNAFILES_H

#include <stdbool.h>
#include "tinycthread.h"

#ifdef __cplusplus
extern "C" {
#endif


/**
 * @brief Get the error number encountered by RnaFile
 * 
 * @return int error number
 */
int rnafgeterrno(void);


/**
 * @brief Error number encountered by RnaFile
 */
#define rnaferrno (rnafgeterrno())


/**
 * @brief Opaque pointer for RnaFile struct
 */
typedef struct RnaFile *RnaFile;


/**
 * @brief Open an fasta/fastq/sequence files fo reading. Mode is used to specify which type
 * of file is used for reading. "a" for fasta, "q" for fastq, "s" for sequence file, and
 * "b" for binary.
 * 
 * @param filename Path to the file you want to open for reading
 * @param mode Type of file being opened
 * @return RnaFile 
 */
RnaFile rnafopen(const char *filename, const char *mode);


/**
 * @brief Close an RnaFile handle.
 * 
 * Removes all allocated resources associated with RnaFile.
 * 
 * @param file RnaFile to close
 */
void rnafclose(RnaFile file);


/**
 * @brief Test if RnaFile is at the end of file.
 * 
 * @param file RnaFile pointer to test
 * @return true if at end of file
 * @return false if not at end of file
 */
bool rnafeof(RnaFile file);


/**
 * @brief Return an allocated string detailing the error encountered from RnaFile
 * 
 * @param rnaferrno_ The rnaferrno variable
 * @return char* String with error description
 * 
 * @note It is necessary to free the string returned from this function as memory is
 * allocated to store the string.
 */
char *rnafstrerror(int _rnaferrno);


/**
 * @brief Read RnaFile sequences into buffer. Will continue reading until it can no
 * longer fit a full sequence within bufsize characters.
 * 
 * @param file    RnaFile pointer to read from
 * @param buffer  Buffer to write sequences to
 * @param bufsize Size of the buffer being passed
 * @return size_t Number of bytes read into buffer. 0 if end of file, or error encountered.
 */
size_t rnafread(RnaFile file, char *buffer, size_t bufsize);


/**
 * @brief Read RnaFile sequences into a buffer. Will continue reading until it can no
 * longer fit a full sequence within butsize characters.
 * 
 * @param file    RnaFile pointer to read from
 * @param buffer  Buffer to write sequences to
 * @param bufsize Size of the buffer being passed
 * @return size_t Number of bytes read into buffer. 0 if end of file, or error encountered.
 * 
 * @note
 * This function does not use a mutex to lock access to the RnaFile buffer. As such, it is not
 * thread safe. Only use in single-threaded applications
 */
size_t rnafread_unlocked(RnaFile file, char *buffer, size_t bufsize);


/** Undocumented read functions. See `rnafread()` for more information. */
size_t rnafqread(RnaFile file, char *buffer, size_t bufsize);
size_t rnafqread_unlocked(RnaFile file, char *buffer, size_t bufsize);
size_t rnafaread(RnaFile file, char *buffer, size_t bufsize);
size_t rnafaread_unlocked(RnaFile file, char *buffer, size_t bufsize);
size_t rnafsread(RnaFile file, char *buffer, size_t bufsize);
size_t rnafsread_unlocked(RnaFile file, char *buffer, size_t bufsize);


/**
 * @brief Read a record's sequence into `buffer`.
 * 
 * `rnafgets()` reads at most `bufsize - 1` characters and stores them into the `buffer` string. Reading
 * stops as soon as the end of the record it is currently processing has been reached or the end of file
 * has been reached. A `'\0'` is appended to `buffer` to ensure the string is null terminated.
 * 
 * @param file    RnaFile to read from
 * @param buffer  Buffer to fill
 * @param bufsize Size of the buffer being passed
 * @return char* Pointer to the record's sequence, or NULL if unable to find the next record.
 */
char *rnafgets(RnaFile file, char *buffer, size_t bufsize);


/**
 * @brief Read a record's sequence into `buffer`.
 * 
 * `rnafgets()` reads at most `bufsize - 1` characters and stores them into the `buffer` string. Reading
 * stops as soon as the end of the record it is currently processing has been reached or the end of file
 * has been reached. A `'\0'` is appended to `buffer` to ensure the string is null terminated.
 * 
 * @param file    RnaFile to read from
 * @param buffer  Buffer to fill
 * @param bufsize Size of the buffer being passed
 * @return char* Pointer to the record's sequence, or NULL if unable to find the next record.
 * 
 * @note
 * This function does not use a mutex to lock access to the RnaFile internal buffer. As such, it is not
 * thread-safe. Only use in single-threaded applications.
 */
char *rnafgets_unlocked(RnaFile file, char *buffer, size_t bufsize);


/** Undocumented gets functions. See rnafgets() for more information. */
char *rnafagets(RnaFile file, char *buffer, size_t bufsize);
char *rnafagets_unlocked(RnaFile file, char *buffer, size_t bufsize);
char *rnafsgets(RnaFile file, char *buffer, size_t bufsize);
char *rnafsgets_unlocked(RnaFile file, char *buffer, size_t bufsize);

#ifdef __cplusplus
}
#endif

#endif
