/*
 Copyright (c) 2017-2018 Tuukka Norri
 This code is licensed under MIT license (see LICENSE for details).
 */

#include <cstdlib>
#include <dispatch/dispatch.h>
#include <founder_sequences/founder_sequences.hh>
#include <founder_sequences/generate_founder_sequences.hh>
#include <iostream>
#include <libbio/assert.hh>
#include <unistd.h>

#ifdef __linux__
#include <pthread_workqueue.h>
#endif

#include "cmdline.h"


namespace lb	= libbio;
namespace fseq	= founder_sequences;


namespace {
	fseq::segment_joining segment_joining_method(enum_segment_joining const sj)
	{
		switch (sj)
		{
			case segment_joining_arg_matching:
				return fseq::segment_joining::MATCHING;
				
			case segment_joining_arg_random:
				return fseq::segment_joining::RANDOM;
				
			case segment_joining__NULL:
			default:
				lb::fail("Unexpected value for structural variant handling.");
				return fseq::segment_joining::MATCHING; // Not reached.
		}
	}
}


int main(int argc, char **argv)
{
	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		exit(EXIT_FAILURE);
	
	if (args_info.segment_length_bound_arg <= 0)
	{
		std::cerr << "Segment length bound must be positive." << std::endl;
		exit(EXIT_FAILURE);
	}

	std::ios_base::sync_with_stdio(false);	// Don't use C style IO after calling cmdline_parser.
	std::cin.tie(nullptr);					// We don't require any input from the user.

#ifndef NDEBUG
	std::cerr << "Assertions have been enabled." << std::endl;
#endif

	// libdispatch on macOS does not need pthread_workqueue.
#ifdef __linux__
	pthread_workqueue_init_np();
#endif

	fseq::generate_founder_sequences(
		args_info.input_arg,
		args_info.segment_length_bound_arg,
		segment_joining_method(args_info.segment_joining_arg),
		args_info.output_segments_arg,
		args_info.single_threaded_flag
	);
		
	cmdline_parser_free(&args_info);
	
	dispatch_main();
	// Not reached b.c. pthread_exit() is eventually called in dispatch_main().
	return EXIT_SUCCESS;
}
