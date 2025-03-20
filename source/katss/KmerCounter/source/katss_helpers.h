#ifndef KATSS_HELPERS_H
#define KATSS_HELPERS_H

#include "katss.h"

/**
 * @brief Parse the KatssOptions options to ensure they are correct.
 * 
 * @param opts Pointer to KatssOptions struct
 * @return int 0 on success, 1 when an error was encountered.
 */
int
katss_parse_options(KatssOptions *opts);


/**
 * @brief Initialize a KatssData struct to the given kmer size
 * 
 * @param kmer Length of k
 * @return KatssData* Pointer to KatssData struct
 */
KatssData *
katss_init_kdata(int kmer);

#endif
