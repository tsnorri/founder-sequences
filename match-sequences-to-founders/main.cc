/*
 Copyright (c) 2018 Tuukka Norri
 This code is licensed under MIT license (see LICENSE for details).
 */

#include <cstdlib>
#include <dispatch/dispatch.h>
#include <founder_sequences/match_founder_sequences.hh>
#include <iostream>
#include <libbio/assert.hh>
#include <libbio/sequence_reader/sequence_reader.hh>
#include <unistd.h>

#ifdef __linux__
#include <pthread_workqueue.h>
#endif

#include "cmdline.h"


namespace fseq	= founder_sequences;
namespace lb	= libbio;
namespace lsr	= libbio::sequence_reader;


namespace {
	
	lsr::input_format input_format(enum_founders_format const fmt)
	{
		switch (fmt)
		{
			case founders_format_arg_FASTA:
				return lsr::input_format::FASTA;
					
			case founders_format_arg_text:
				return lsr::input_format::TEXT;
				
			case founders_format_arg_listMINUS_file:
				return lsr::input_format::LIST_FILE;
			
			case founders_format__NULL:
			default:
				libbio_fail("Unexpected input format.");
				return lsr::input_format::LIST_FILE;
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

	fseq::match_founder_sequences(
		args_info.sequences_arg,
		args_info.founders_arg,
		input_format(args_info.founders_format_arg),
		args_info.single_threaded_flag
	);
	
	cmdline_parser_free(&args_info);
	
	dispatch_main();
	// Not reached b.c. pthread_exit() is eventually called in dispatch_main().
	return EXIT_SUCCESS;
}
