/*
 Copyright (c) 2017-2018 Tuukka Norri
 This code is licensed under MIT license (see LICENSE for details).
 */

#include <cstdlib>
#include <dispatch/dispatch.h>
#include <founder_sequences/founder_sequences.hh>
#include <founder_sequences/generate_context.hh>
#include <iostream>
#include <libbio/assert.hh>
#include <unistd.h>

#ifdef __linux__
#include <pthread_workqueue.h>
#endif

#include "cmdline.h"


namespace lb	= libbio;
namespace lsr	= libbio::sequence_reader;
namespace fseq	= founder_sequences;


namespace {
	
	fseq::segment_joining segment_joining_method(enum_segment_joining const sj)
	{
		switch (sj)
		{
			case segment_joining_arg_bipartiteMINUS_matching:
				return fseq::segment_joining::BIPARTITE_MATCHING;

			case segment_joining_arg_greedy:
				return fseq::segment_joining::GREEDY;

			case segment_joining_arg_random:
				return fseq::segment_joining::RANDOM;

			case segment_joining__NULL:
			default:
				libbio_fail("Unexpected value for structural variant handling.");
				return fseq::segment_joining::GREEDY; // Not reached.
		}
	}
	
	
	lsr::input_format input_file_format(enum_input_format const fmt)
	{
		switch (fmt)
		{
			case input_format_arg_FASTA:
				return lsr::input_format::FASTA;
			
			case input_format_arg_listMINUS_file:
				return lsr::input_format::LIST_FILE;
			
			case input_format__NULL:
			default:
				libbio_fail("Unexpected value for input_format");
				return lsr::input_format::LIST_FILE; // Not reached.
		}
	}
}


int main(int argc, char **argv)
{
	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		exit(EXIT_FAILURE);

	std::ios_base::sync_with_stdio(false);	// Don't use C style IO after calling cmdline_parser.

#ifndef NDEBUG
	std::cerr << "Assertions have been enabled." << std::endl;
#endif

	// Handle the command line arguments.
	if (args_info.print_invocation_flag)
	{
		std::cerr << "Invocation:";
		for (std::size_t i(0); i < argc; ++i)
			std::cerr << ' ' << argv[i];
		std::cerr << std::endl;
	}
	
	auto const mode(fseq::running_mode::GENERATE_FOUNDERS);
	
	if (args_info.segment_length_bound_given)
	{
		if (args_info.segment_length_bound_arg <= 0)
		{
			std::cerr << "Segment length bound must be positive." << std::endl;
			exit(EXIT_FAILURE);
		}
	}
	else
	{
		std::cerr << "Segment length bound needs to be specified when generating a segmentation." << std::endl;
		exit(EXIT_FAILURE);
	}
	
	if (! (0 <= args_info.random_seed_arg && args_info.random_seed_arg <= std::numeric_limits <std::uint_fast32_t>::max()))
	{
		std::cerr << "Random seed out of bounds." << std::endl;
		exit(EXIT_FAILURE);
	}
	
	if (args_info.pbwt_sample_rate_arg <= 0)
	{
		std::cerr << "PBWT sample rate multiplier must be non-negative." << std::endl;
		exit(EXIT_FAILURE);
	}
	
	auto const segment_joining(segment_joining_method(args_info.segment_joining_arg));
	
	// Instantiate the controller class and run.
	{
		// Deallocates itself with a callback.
		auto *ctx(new fseq::generate_context(
			mode,
			args_info.segment_length_bound_arg,
			segment_joining,
			fseq::bipartite_set_scoring::INTERSECTION,
			args_info.pbwt_sample_rate_arg,
			args_info.random_seed_arg,
			args_info.single_threaded_flag
		));
		
		ctx->prepare(
			nullptr,
			nullptr,
			args_info.output_founders_arg,
			args_info.output_segments_arg
		);
		ctx->load_and_generate(
			args_info.input_arg,
			input_file_format(args_info.input_format_arg)
		);
	}
	
	// Everything in args_info should have been copied by now, so it is no longer needed.
	cmdline_parser_free(&args_info);
	
	dispatch_main();
	// Not reached b.c. pthread_exit() is eventually called in dispatch_main().
	return EXIT_SUCCESS;
}
