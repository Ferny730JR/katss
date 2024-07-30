#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "enrichments.h"
#include "hash_functions.h"
#include "ikke_cmdl.h"
#include "memory_utils.h"
#include "string_utils.h"

typedef struct Options {
	char *test_file;
	char *ctrl_file;
	FILE *out_file;

	unsigned int kmer;
	int iterations;
	char file_delimiter;
	bool no_log;

	bool enrichments;
	bool probabilistic;
} Options;

char 
delimiter_to_char(char *user_delimiter);

int
set_out_file(Options *opt, char *name);

void
enrichments_to_file(KatssEnrichments *enrichments, Options *opt);

void
options_free(Options *opt);

void
init_default_options(Options *opt)
{
	opt->test_file = NULL;
	opt->ctrl_file = NULL;
	opt->out_file  = NULL;

	opt->kmer = 5;
	opt->iterations = 1;
	opt->file_delimiter = ',';
	opt->no_log = false;

	opt->enrichments = false;
	opt->probabilistic = false;
}


int
main(int argc, char *argv[])
{
	struct ikke_args_info    args_info;
	Options                  opt;

	init_default_options(&opt);

	/* Parse command line options */
	if(ikke_cmdline_parser(argc, argv, &args_info) != 0)
		goto exit_error;
	
	if(!args_info.test_given && !args_info.control_given) {
		printf("ikke -c [control.fastq.gz] -t [test.fastq.gz] [OPTIONS]\n");
		printf("Try 'ikke --help' for more information.\n");
		goto cleanup_args;
	}

	if(args_info.test_given)
		opt.test_file = strdup(args_info.test_arg);
	else {
		error_message("You need to provide a 'test' file.");
		goto cleanup_args;
	}

	if(args_info.control_given)
		opt.ctrl_file = strdup(args_info.control_arg);

	if(args_info.file_delimiter_given)
		opt.file_delimiter = delimiter_to_char(args_info.file_delimiter_arg);

	if(set_out_file(&opt, args_info.output_arg) != 0)
		goto cleanup_args;

	opt.kmer          = args_info.kmer_arg;
	opt.iterations    = args_info.iterations_arg;
	opt.no_log        = (bool)args_info.no_log_flag;
	opt.enrichments   = (bool)args_info.enrichments_flag;
	opt.probabilistic = (bool)args_info.independent_probs_flag;

	/* Make sure options were set correctly */
	if(opt.kmer > 16) {
		error_message("Option 'kmer' must be a value between 0 and 16\n"
		              "Actual: %u", opt.kmer);
		goto cleanup_args;
	}

	if(opt.ctrl_file == NULL && !opt.probabilistic) {
		error_message("You need to provide a control file");
		goto cleanup_args;
	}

	if(opt.probabilistic && opt.ctrl_file != NULL)
		warning_message("Ignoring control file: %s", opt.ctrl_file);

	if(opt.iterations < 1) {
		warning_message("Iterations can not be below 1. Defaulting to 1.");
		opt.iterations = 1;
	}

	if(opt.iterations > 1ULL << (2*opt.kmer)) {
		warning_message("Iterations can not be larger than 4^kmer. Deafulting to 1.");
		opt.iterations = 1;
	}

	ikke_cmdline_parser_free(&args_info);

	/* Begin actual computations */
	KatssEnrichments *enr;
	if(opt.enrichments && opt.probabilistic)
		enr = katss_prob_enrichments(opt.test_file, opt.kmer, !opt.no_log);
	else if(opt.enrichments)
		enr = katss_enrichments(opt.test_file, opt.ctrl_file, opt.kmer, !opt.no_log);
	else if(opt.probabilistic)
		enr = katss_prob_ikke(opt.test_file, opt.kmer, opt.iterations, !opt.no_log);
	else
		enr = katss_ikke(opt.test_file, opt.ctrl_file, opt.kmer, opt.iterations, !opt.no_log);

	/* Check results */
	if(enr == NULL)
		goto cleanup_opts;
	
	enrichments_to_file(enr, &opt);

	/* Everything seemed to work! Cleanup and return */
	katss_free_enrichments(enr);
	options_free(&opt);
	return 0;

cleanup_args:
	ikke_cmdline_parser_free(&args_info);
cleanup_opts:
	options_free(&opt);
exit_error:
	exit(EXIT_FAILURE);
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
		name = "motif";

	char *out_filename;
	switch(opt->file_delimiter) {
		case ',':  out_filename = concat(name, ".csv"); break;
		case '\t': out_filename = concat(name, ".tsv"); break;
		default:   out_filename = concat(name, ".dsv"); break;
	}

	opt->out_file = fopen(out_filename, "w");
	free(out_filename);

	return opt->out_file == NULL;
}

void
enrichments_to_file(KatssEnrichments *enrichments, Options *opt)
{
	fprintf(opt->out_file, "kmer%crval\n", opt->file_delimiter);
	char kseq[17];
	for(uint64_t i = 0; i < enrichments->num_enrichments; i++) {
		double rval = enrichments->enrichments[i].enrichment;
		if(isnan(rval))
			continue;
		katss_unhash(kseq, enrichments->enrichments[i].key, opt->kmer, true);
		fprintf(opt->out_file, "%s%c%f\n", kseq, opt->file_delimiter, rval);
	}
}

void
options_free(Options *opt)
{
	if(opt->test_file)
		free(opt->test_file);
	if(opt->ctrl_file)
		free(opt->ctrl_file);
	if(opt->out_file)
		fclose(opt->out_file);
}
