/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef FOUNDER_SEQUENCES_GENERATE_FOUNDER_SEQUENCES_HH
#define FOUNDER_SEQUENCES_GENERATE_FOUNDER_SEQUENCES_HH


namespace founder_sequences {
	
	void generate_founder_sequences(
		char const *input_path,
		input_format const input_file_format,
		std::size_t const segment_length,
		segment_joining const segment_joining_method,
		char const *output_segments_path,
		char const *output_founders_path,
		bool const use_single_thread
	);
}

#endif
