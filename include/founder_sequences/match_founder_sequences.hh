/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef FOUNDER_SEQUENCES_MATCH_FOUNDER_SEQUENCES_HH
#define FOUNDER_SEQUENCES_MATCH_FOUNDER_SEQUENCES_HH


namespace founder_sequences {
	
	void match_founder_sequences(
		char const *sequences_list_path,
		char const *founders_list_path,
		bool const single_threaded
	);
}

#endif
