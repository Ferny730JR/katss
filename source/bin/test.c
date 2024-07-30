#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

#include "enrichments.h"

int
main(int argc, char *argv[])
{
    if(argc != 3) {
        printf("You need to pass a filepath and kmer, in that order!\n");
        return 1;
    }

    const char *file = argv[1];

    unsigned int kmer = atoi(argv[2]);
    KatssCounter *counter = katss_count_kmers(file, kmer);
    if(counter == NULL) {
        printf("Failed to get counts!");
        return 1;
    }
    printf("Got counts!\n");
    katss_free_counter(counter);

    return 0;
}
