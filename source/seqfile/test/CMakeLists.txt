add_executable(test-seqf tests.c)
target_link_libraries(test-seqf seqf_static)
target_compile_definitions(test-seqf PRIVATE
    EXAMPLE_FASTA=${CMAKE_CURRENT_SOURCE_DIR}/example_files/example.fasta
    EXAMPLE_FASTA_GZ=${CMAKE_CURRENT_SOURCE_DIR}/example_files/example.fasta.gz
    EXAMPLE_FASTQ=${CMAKE_CURRENT_SOURCE_DIR}/example_files/example.fastq
    EXAMPLE_FASTQ_GZ=${CMAKE_CURRENT_SOURCE_DIR}/example_files/example.fastq.gz
    EXAMPLE_READS=${CMAKE_CURRENT_SOURCE_DIR}/example_files/example.reads
    EXAMPLE_READS_GZ=${CMAKE_CURRENT_SOURCE_DIR}/example_files/example.reads.gz
)
