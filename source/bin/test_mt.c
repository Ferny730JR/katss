#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

#include "enrichments.h"

int
main(int argc, char *argv[])
{
    if(argc != 4) {
        printf("You need to pass a filepath, kmer, and threads in that order!\n");
        return 1;
    }

    const char *file = argv[1];

    unsigned int kmer = atoi(argv[2]);
	int threads = atoi(argv[3]);

    KatssCounter *counter = katss_count_kmers_mt(file, kmer, threads);
    if(counter == NULL) {
        printf("Failed to get counts!");
        return 1;
    }
    printf("Got counts!\n");
    katss_free_counter(counter);

    return 0;
}
