#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "enrichments.h"
#include "bootstrap.h"
#include "hash_functions.h"
#include "ikke_cmdl.h"
#include "memory_utils.h"
#include "string_utils.h"

typedef struct Options {
	char *test_file;    /** Name of test file */
	char *ctrl_file;    /** Name of control file */
	FILE *out_file;     /** Output file to write to */

	int  kmer;          /** Length of k-mer to count */
	int  iterations;    /** Number of ikke iterations */
	int  threads;       /** Number of threads to use */
	char delimiter;     /** File delimiter for output file */
	bool no_log;        /** Don't normalize outputs to log2 */

	bool enrichments;   /** Compute regular enrichments */
	bool probabilistic; /** Enable probabilistic enrichments */
	bool bootstrap;     /** Enable bootstrap */
	int  bs_runs;       /** Bootstrap iterations to perform */
	int  sample;        /** Percent of file to sample */
} Options;

void
process_bootstrap(Options *opts);

char 
delimiter_to_char(char *user_delimiter);

int
set_out_file(Options *opt, char *name);

void
enrichments_to_file(KatssEnrichments *enrichments, Options *opt);

void
bootstrap_to_file(KatssBootstrap *bootstrap, Options *opts);

void
options_free(Options *opt);

void
init_default_options(Options *opt)
{
	opt->test_file = NULL;
	opt->ctrl_file = NULL;
	opt->out_file  = NULL;

	opt->kmer       = 5;
	opt->iterations = 1;
	opt->threads    = 1;
	opt->delimiter  = ',';
	opt->no_log     = false;

	opt->enrichments   = false;
	opt->probabilistic = false;
	opt->bootstrap     = false;
	opt->bs_runs       = 10;
	opt->sample        = 10;
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

	if(args_info.delimiter_given)
		opt.delimiter = delimiter_to_char(args_info.delimiter_arg);

	if(set_out_file(&opt, args_info.output_arg) != 0)
		goto cleanup_args;

	opt.kmer          = args_info.kmer_arg;
	opt.iterations    = args_info.iterations_arg;
	opt.threads       = args_info.threads_arg;
	opt.sample        = args_info.sample_arg;
	opt.bs_runs       = args_info.bootstrap_arg;
	opt.bootstrap     = args_info.bootstrap_given;
	opt.no_log        = (bool)args_info.no_log_flag;
	opt.enrichments   = (bool)args_info.enrichments_flag;
	opt.probabilistic = (bool)args_info.independent_probs_flag;

	/* Make sure options were set correctly */
	if(opt.kmer < 1 || 16 < opt.kmer) {
		error_message("Option 'kmer' must be a value between 1 and 16. "
		              "Given: %d", opt.kmer);
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

	if(!opt.bootstrap && args_info.sample_given)
		warning_message("Ignoring sample. Bootstrap must be enable to sub-sample file.");

	if(opt.bootstrap && (opt.sample < 1 || 100 < opt.sample)) {
		warning_message("Invalid sample: cannot sample %d%% of the file. "
		                "Sample must be a number between 1 and 100. Defaulting "
						"to 10.");
		opt.sample = 10;
	}

	ikke_cmdline_parser_free(&args_info);

	/* Begin actual computations */
	if(opt.bootstrap) {
		process_bootstrap(&opt);
	} else {
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
		katss_free_enrichments(enr);
	}

	/* Everything seemed to work! Cleanup and return */
	options_free(&opt);
	return 0;

cleanup_args:
	ikke_cmdline_parser_free(&args_info);
cleanup_opts:
	options_free(&opt);
exit_error:
	exit(EXIT_FAILURE);
}

void
process_bootstrap(Options *opts)
{
	/* Sanity check... */
	if(!opts->bootstrap) {
		error_message("Bootstrap not enabled, even though it appeared that it was.");
		return;
	}

	/* Declare bootstrap and its options */
	KatssBootstrap *bootstrap;
	KatssOptions bootstrap_opts;

	/* Initialize default options, then populate it with params */
	katss_init_default_opts(&bootstrap_opts);
	bootstrap_opts.algo          = opts->enrichments ? enrichments : ikke;
	bootstrap_opts.bs_iters      = opts->bs_runs;
	bootstrap_opts.ikke_iters    = opts->iterations;
	bootstrap_opts.kmer          = opts->kmer;
	bootstrap_opts.probabilistic = opts->probabilistic;
	bootstrap_opts.sample        = opts->sample;
	bootstrap_opts.threads       = opts->threads;
	
	/* Process bootstrap */
	bootstrap = katss_bootstrap(opts->test_file, opts->ctrl_file, &bootstrap_opts);
	if(bootstrap == NULL)
		return; /* error encountered, don't output to file */
	
	bootstrap_to_file(bootstrap, opts);
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
	switch(opt->delimiter) {
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
	fprintf(opt->out_file, "kmer%crval\n", opt->delimiter);
	char kseq[17];
	for(uint64_t i = 0; i < enrichments->num_enrichments; i++) {
		double rval = enrichments->enrichments[i].enrichment;
		if(isnan(rval))
			continue;
		katss_unhash(kseq, enrichments->enrichments[i].key, opt->kmer, true);
		fprintf(opt->out_file, "%s%c%f\n", kseq, opt->delimiter, rval);
	}
}

void
bootstrap_to_file(KatssBootstrap *bootstrap, Options *opts)
{
	fprintf(opts->out_file, "kmer%cmean%cstdev\n", opts->delimiter, opts->delimiter);
	char kseq[17];
	for(uint64_t i = 0; i < bootstrap->total; i++) {
		katss_unhash(kseq, bootstrap->data[i].kmer_hash, opts->kmer, true);
		fprintf(opts->out_file, "%s%c%f%c%f\n", kseq, opts->delimiter,
		 bootstrap->data[i].mean, opts->delimiter, bootstrap->data[i].stdev);
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
