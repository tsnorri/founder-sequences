/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef FOUNDER_SEQUENCES_SEGMENTATION_CONTEXT_HH
#define FOUNDER_SEQUENCES_SEGMENTATION_CONTEXT_HH

#include <founder_sequences/founder_sequences.hh>


namespace founder_sequences {
	
	struct segmentation_context
	{
		virtual ~segmentation_context() {}
		virtual std::uint32_t max_segment_size() const = 0;
	};
	
	
	struct segmentation_context_delegate
	{
		virtual ~segmentation_context_delegate() {}
		virtual alphabet_type const &alphabet() const = 0;
		virtual sequence_vector const &sequences() const = 0;
		virtual bipartite_set_scoring bipartite_set_scoring_method() const = 0;
		virtual bool should_run_single_threaded() const = 0;
		
		virtual std::uint32_t sequence_count() const = 0;
		
		virtual std::ostream &sequence_output_stream() = 0;
		virtual std::ostream &segments_output_stream() = 0;
	};
}

#endif
