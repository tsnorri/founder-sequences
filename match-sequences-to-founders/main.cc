/*
 Copyright (c) 2018 Tuukka Norri
 This code is licensed under MIT license (see LICENSE for details).
 */

#include <cstdlib>
#include <dispatch/dispatch.h>
#include <founder_sequences/match_founder_sequences.hh>
#include <iostream>
#include <libbio/assert.hh>
#include <unistd.h>

#ifdef __linux__
#include <pthread_workqueue.h>
#endif

#include "cmdline.h"


namespace lb	= libbio;
namespace fseq	= founder_sequences;


int main(int argc, char **argv)
{
	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		exit(EXIT_FAILURE);
	
	std::ios_base::sync_with_stdio(false);	// Don't use C style IO after calling cmdline_parser.
	std::cin.tie(nullptr);					// We don't require any input from the user.

#ifndef NDEBUG
	std::cerr << "Assertions have been enabled." << std::endl;
#endif

	// libdispatch on macOS does not need pthread_workqueue.
#ifdef __linux__
	pthread_workqueue_init_np();
#endif
	
	fseq::match_founder_sequences(
		args_info.sequences_list_arg,
		args_info.founders_list_arg,
		args_info.single_threaded_flag
	);
	
	cmdline_parser_free(&args_info);
	
	dispatch_main();
	// Not reached b.c. pthread_exit() is eventually called in dispatch_main().
	return EXIT_SUCCESS;
}
