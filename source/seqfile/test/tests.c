#include <stdio.h>

#include "minunit.h"

#include "seqf_core.h"

#define STRINGIZE(arg) #arg
#define TXT2STR(arg) STRINGIZE(arg)

static UTEST_TYPE
test_seqfopen(void)
{
	init_unit_tests("Testing seqfopen");
	SeqFile file;

#define assert_file(file, exp_compression, exp_type) \
	file != NULL && \
	((seqf_statep)file)->compression == exp_compression && \
	((seqf_statep)file)->type == exp_type && \
	((seqf_statep)file)->mutex_is_init

	file = seqfopen(TXT2STR(EXAMPLE_FASTA), "a");
	mu_assert("Open fasta file", assert_file(file, PLAIN, 'a'));
	seqfclose(file);

	file = seqfopen(TXT2STR(EXAMPLE_FASTA_GZ), "a");
	mu_assert("Open compressed fasta file", assert_file(file, GZIP, 'a'));
	seqfclose(file);

	file = seqfopen(TXT2STR(EXAMPLE_FASTQ), "q");
	mu_assert("Open fastq file", assert_file(file, PLAIN, 'q'));
	seqfclose(file);

	file = seqfopen(TXT2STR(EXAMPLE_FASTQ_GZ), "q");
	mu_assert("Open compressed fastq file", assert_file(file, GZIP, 'q'));
	seqfclose(file);

	file = seqfopen(TXT2STR(EXAMPLE_READS), "s");
	mu_assert("Open sequences file", assert_file(file, PLAIN, 's'));
	seqfclose(file);

	file = seqfopen(TXT2STR(EXAMPLE_READS_GZ), "s");
	mu_assert("Open compressed sequences file", assert_file(file, GZIP, 's'));
	seqfclose(file);

	/* Check that return NULL on files that don't exist */
	file = seqfopen("non-existent-dir/non-existent-file", NULL);
	mu_assert("Opening non-existent file", file == NULL);
	seqfclose(file);

	file = seqfopen("example_files/example.reads", "wrong!");
	mu_assert("Open file with unsupported mode", file == NULL);

#undef assert_file
	unit_tests_end;
}

static UTEST_TYPE
test_seqfclose(void)
{
	init_unit_tests("Testing seqfclose");
	mu_assert("Close null file", seqfclose(NULL) == 1);
	mu_assert("Close opened file", seqfclose(seqfopen(TXT2STR(EXAMPLE_READS),NULL))==0);

	unit_tests_end;
}

static UTEST_TYPE
test_seqferrno(void)
{
	init_unit_tests("Testing seqferrno");

	mu_assert("seqferrno has not been set", seqferrno == 0);

	seqfopen("non-existent-path/no-existent-file", NULL);
	mu_assert("Opening file that does not exist", seqferrno == 1);

	seqfopen(TXT2STR(EXAMPLE_READS), "wrong mode!");
	mu_assert("Error code for wrong mode", seqferrno == 3);

	unit_tests_end;
}

static UTEST_TYPE
test_seqfgetc(void)
{
	init_unit_tests("Testing seqfgetc");

	/* Prepare variables to use */
	bool passed = false;
	SeqFile file = seqfopen(TXT2STR(EXAMPLE_READS), NULL);

	passed = seqfgetc(file) == 'G';
	mu_assert("getc call on uncompressed file", passed);

	seqfrewind(file); passed = true;
	char *reads_seq = "GCATACGGTGAAAGCTCAGCTTTCCAGCGCTGCTTTACAGTTGGCACGATTAACCCAAGAACGTTATTTCTGTCAAATTTTGAGTTGGTTGTGGGCAAGG";
	char *tmp = reads_seq;
	for(int i=0; i<100; i++) {
		if(*tmp++ != seqfgetc(file)) {
			passed = false;
			break;
		}
	}
	mu_assert("Many seqfgetc calls", passed);

	seqfclose(file);
	file = seqfopen(TXT2STR(EXAMPLE_READS_GZ), NULL);
	passed = seqfgetc(file) == 'G';
	mu_assert("getc call on compressed file", passed);

	seqfrewind(file);
	passed = true;
	tmp = reads_seq;
	for(int i=0; i<100; i++) {
		if(*tmp++ != seqfgetc(file)) {
			passed = false;
			break;
		}
	}
	seqfclose(file);
	mu_assert("Many seqfgetc calls on compressed file", passed);

	char buf[300]; tmp=buf;
	FILE *fp = fopen(TXT2STR(EXAMPLE_FASTA), "r");
	file = seqfopen(TXT2STR(EXAMPLE_FASTA), NULL);
	int nread = fread(buf, 1, 290, fp);
	for(int i=0; i<nread; i++) {
		if(*tmp++ != seqfgetc(file)) {
			passed = false;
			break;
		}
	}
	fclose(fp);
	mu_assert("seqfgetc entire file", passed);

	passed = seqfgetc(file) == EOF;
	mu_assert("seqfgetc when reached end of file", passed);
	seqfclose(file);

	unit_tests_end;
}

static void all_tests() {
	init_run_test;

	/* Begin tests */
	mu_run_test(test_seqfopen);
	mu_run_test(test_seqfclose);
	mu_run_test(test_seqferrno);
	mu_run_test(test_seqfgetc);

	/* End of tests */
	run_test_end;
}

int main(void) {
	all_tests();
	return 0;
}
