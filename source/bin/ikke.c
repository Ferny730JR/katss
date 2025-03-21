#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "hash_functions.h"
#include "ikke_cmdl.h"
#include "memory_utils.h"
#include "string_utils.h"

#include "katss.h"

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
	bool shuffle;       /** Shuffle the sequences */
	int  klet;          /** Length of k-let to preserve during shuffling */
	bool probabilistic; /** Enable probabilistic enrichments */
	bool bootstrap;     /** Enable bootstrap */
	int  bs_runs;       /** Bootstrap iterations to perform */
	int  sample;        /** Percent of file to sample */
} Options;

char 
delimiter_to_char(char *user_delimiter);

int
set_out_file(Options *opt, char *name);

void
katssdata_to_file(KatssData *data, Options *opt);

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
	opt.klet          = args_info.klet_arg;
	opt.shuffle       = (bool)args_info.shuffle_flag;
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

	if(opt.no_log)
		warning_message("ikke: option --no-log is being ignored. Values are no longer normalized to log2");
	opt.no_log = true;

	ikke_cmdline_parser_free(&args_info);

	/* Specify all the katss options */
	KatssOptions katss_opts;
	katss_init_options(&katss_opts);

	katss_opts.kmer = opt.kmer;
	katss_opts.iters = opt.iterations;
	katss_opts.threads = opt.threads;
	katss_opts.normalize = !opt.no_log;
	katss_opts.sort_enrichments = true;
	katss_opts.bootstrap_iters = opt.bootstrap ? opt.bs_runs : 0;
	katss_opts.bootstrap_sample = opt.sample*1000;
	katss_opts.probs_ntprec = opt.klet;
	katss_opts.seed = 1;
	if(opt.probabilistic && opt.shuffle) {
		katss_opts.probs_algo = KATSS_PROBS_BOTH;
	} else if(opt.probabilistic) {
		katss_opts.probs_algo = KATSS_PROBS_REGULAR;
	} else if(opt.shuffle) {
		katss_opts.probs_algo = KATSS_PROBS_USHUFFLE;
	} else {
		katss_opts.probs_algo = KATSS_PROBS_NONE;
	}

	/* Begin actual computations */
	KatssData *data = NULL;
	if(opt.enrichments) {
		data = katss_enrichment(opt.test_file, opt.ctrl_file, &katss_opts);
	} else {
		data = katss_ikke(opt.test_file, opt.ctrl_file, &katss_opts);
	}

	/* Failed to get enrichments */
	if(data == NULL)
		goto cleanup_opts;

	/* Output data into a file */
	katssdata_to_file(data, &opt);

	/* Everything seemed to work! Cleanup and return */
	katss_free_kdata(data);
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
katssdata_to_file(KatssData *data, Options *opt)
{
	/* Print the header */
	if(opt->bootstrap) {
		fprintf(opt->out_file, "kmer%crval%cstdev%cpval\n",
		  opt->delimiter, opt->delimiter, opt->delimiter);
	} else {
		fprintf(opt->out_file, "kmer%crval\n",opt->delimiter);
	}

	char kseq[32];
	for(uint64_t i=0; i<data->num_kmers; i++) {
		double rval = data->kmers[i].rval;
		if(isnan(rval))
			continue;
		katss_unhash(kseq, data->kmers[i].kmer, opt->kmer, true);
		if(opt->bootstrap) {
			fprintf(opt->out_file, "%s%c%f%c%f%c%E\n", kseq, opt->delimiter,
			  rval, opt->delimiter, data->kmers[i].stdev, opt->delimiter,
			  data->kmers[i].pval);
		} else {
			fprintf(opt->out_file, "%s%c%f\n", kseq, opt->delimiter, rval);
		}
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
