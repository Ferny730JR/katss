#ifndef KATSS_SEQSEQ
#define KATSS_SEQSEQ

/**
 * @brief Find a sequence in a sequence.
 * 
 * @param seq Sequence to search
 * @param pat Pattern to search for
 * @return char* First occurrence of `pat` in `seq`, otherwise NULL
 */
char *seqseq(const char *seq, const char *pat);


/**
 * @brief Find a sequence in a sequence. Returns the beginning of the line that
 * `pat` was found within `seq`.
 * 
 * @param seq Sequence to search
 * @param pat Pattern to search for
 * @return char* Beginning of line that `pat` was found in. NULL if not found.
 */
char *seqlseq(const char *seq, const char *pat);


/**
 * @brief Find a sequence in a fasta formatted sequence. Doesn't search for `pat`
 * in lines in `seq` that begin with `'>'`.
 * 
 * @param seq Fasta formatted sequence to search
 * @param pat Pattern to search for
 * @return char* First occurrence of `pat` in `seq`, otherwise NULL
 */
char *seqseqa(const char *seq, const char *pat);


/**
 * @brief Find the Beginning of a sequence in a fasta formatted sequence.
 * 
 * @param seq Fasta formatted sequence to search
 * @param pat Pattern to search for
 * @return char* Beginning of sequence in which `pat` was found. NULL if not found.
 */
char *seqlseqa(const char *seq, const char *pat);


/**
 * @brief Find a sequence in a fastq formatted sequence. 
 * 
 * Doesn't search for `pat` in lines in `seq` that begin with  `'\@'`,
 * or lines that begin with `'+'`, including the following line.
 * 
 * @param seq Fastq formatted sequence to search
 * @param pat Pattern to search for
 * @return char* First occurrence of `pat` in `seq`, otherwise NULL
 */
char *seqseqq(const char *seq, const char *pat);


/**
 * @brief Find the line of a sequence in a fastq formatted sequence.
 * 
 * @param seq Fastq formatted sequence to search
 * @param pat Pattern to search for
 * @return char* Beginning of line where `pat` was found. NULL if not found.
 */
char *seqlseqq(const char *seq, const char *pat);

#endif // KATSS_SEQSEQ
