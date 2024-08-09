#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "kstruct_cmdl.h"

#include "structure.h"
#include "string_utils.h"
#include "memory_utils.h"

typedef struct Options {
	char *test_file;   /** Name of test file */
	char *ctrl_file;   /** Name of control file */
	char *out_file;    /** Output file to write to */

	int  kmer;         /** Length of k-mer to count */
	int  threads;      /** Number of threads to use */
	char delimiter;    /** File delimiter for output file */

	bool seq_windows;  /** Use the sliding window algorithm? */
	int  window_size;  /** Length of window to obtain from sequence */
} Options;

char 
delimiter_to_char(char *user_delimiter);

int
set_out_file(Options *opt, char *name);

void
free_options(Options *opt);

void
init_default_options(Options *opts)
{
	opts->test_file = NULL;
	opts->ctrl_file = NULL;
	opts->out_file  = NULL;

	opts->kmer      = 5;
	opts->threads   = 1;
	opts->delimiter = ',';

	opts->seq_windows = 20;
}

int
main(int argc, char *argv[])
{
	struct kstruct_args_info args_info;
	Options                  opts;

	init_default_options(&opts);

	/* Parse command line options */
	if(kstruct_cmdline_parser(argc, argv, &args_info) != 0)
		goto exit_failure;

	if(!args_info.test_arg || !args_info.control_arg) {
		printf("kstruct -t [test.fastq.gz] -c [ctrl.fastq.gz] [OPTIONS]\n");
		printf("Try 'kstruct --help' for more information.\n");
		goto cleanup_args;
	}

	if(args_info.test_given)
		opts.test_file = strdup(args_info.test_arg);
	else {
		error_message("You need to provide a 'test' file.");
		goto cleanup_args;
	}

	if(args_info.control_given)
		opts.ctrl_file = strdup(args_info.control_arg);
	else {
		error_message("You need to provide a 'control' file.");
		goto cleanup_args;
	}

	if(args_info.delimiter_given)
		opts.delimiter = delimiter_to_char(args_info.delimiter_arg);

	if(set_out_file(&opts, args_info.output_arg) != 0)
		goto cleanup_args;

	opts.kmer          = args_info.kmer_arg;
	opts.threads       = args_info.threads_arg;
	opts.seq_windows   = args_info.seq_windows_given;
	opts.window_size   = args_info.seq_windows_arg;

	/* Make sure options were set correctly */
	if(opts.kmer < 1 || 16 < opts.kmer) {
		error_message("Option 'kmer' must be a value between 1 and 16. "
		              "Given: %d", opts.kmer);
		goto cleanup_args;
	}

	if(opts.window_size < 1) {
		error_message("Option 'seq-windows' must be at least one. "
		              "Given: %d", opts.window_size);
	}

	kstruct_cmdline_parser_free(&args_info);

	/* Setup options for katss_bpp */
	BppOptions bpp_opts;
	katss_bpp_init_default_opts(&bpp_opts);
	bpp_opts.kmer        = opts.kmer;
	bpp_opts.seq_windows = opts.seq_windows;
	bpp_opts.window_size = opts.window_size;
	bpp_opts.threads     = opts.threads;

	/* Calculate structural preference */
	kmerHashTable *enrichments;
	enrichments = katss_bpp(opts.test_file, opts.ctrl_file, &bpp_opts);
	if(enrichments == NULL)
		goto cleanup_opts;
	
	/* Output to file */
	kmerHashTable_to_file(enrichments, opts.out_file, opts.delimiter);

	/* Cleanup and return */
	katss_free_bpp(enrichments);
	free_options(&opts);
	return 0;
cleanup_args:
	kstruct_cmdline_parser_free(&args_info);
cleanup_opts:
	free_options(&opts);
exit_failure:
	return 1;
}

char 
delimiter_to_char(char *user_delimiter)
{
	char delimiter;
	char provided_val = user_delimiter[0];
	if(strlen(user_delimiter)>1) {
		provided_val = 0;
	}

	switch(provided_val) {
		case ',':
			delimiter = ',';
			break;
		case 't':
			delimiter = '\t';
			break;
		case ':':
			delimiter = ':';
			break;
		case '|':
			delimiter = '|';
			break;
		case ' ':
			delimiter = ' ';
			break;
		default:
			warning_message("Provided delimiter '%s' is not valid. Defaulting "
			 "to ','", user_delimiter);
			delimiter = ',';
	}

	return delimiter;
}

int
set_out_file(Options *opt, char *name)
{
	if(name == NULL)
		name = "rna";

	char *out_filename;
	switch(opt->delimiter) {
		case ',':  out_filename = concat(name, ".csv"); break;
		case '\t': out_filename = concat(name, ".tsv"); break;
		default:   out_filename = concat(name, ".dsv"); break;
	}

	opt->out_file = out_filename;
	return opt->out_file == NULL;
}

void
free_options(Options *opt)
{
	if(opt->test_file)
		free(opt->test_file);
	if(opt->ctrl_file)
		free(opt->ctrl_file);
	if(opt->out_file)
		free(opt->out_file);
}
