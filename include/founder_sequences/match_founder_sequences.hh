/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef FOUNDER_SEQUENCES_MATCH_FOUNDER_SEQUENCES_HH
#define FOUNDER_SEQUENCES_MATCH_FOUNDER_SEQUENCES_HH

#include <libbio/sequence_reader/sequence_reader.hh>


namespace founder_sequences {
	
	void match_founder_sequences(
		char const *sequences_list_path,
		char const *founders_path,
		libbio::sequence_reader::input_format const founders_format,
		std::size_t const min_segment_length,
		bool const single_threaded
	);
}

#endif
