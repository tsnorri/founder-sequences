/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef FOUNDER_SEQUENCES_MATCHER_HH
#define FOUNDER_SEQUENCES_MATCHER_HH

#include <founder_sequences/founder_sequences.hh>
#include <founder_sequences/segmentation_dp_arg.hh>


namespace founder_sequences {
	
	struct matcher
	{
		virtual ~matcher() {}
		virtual void match() = 0;
		virtual void output_segments(std::ostream &stream, sequence_vector const &sequences) = 0;
	};


	struct matcher_delegate
	{
		virtual ~matcher_delegate() {}
	
		virtual libbio::dispatch_ptr <dispatch_queue_t> producer_queue() const = 0;	// Copy b.c. the pointer manages a reference counted object.
		virtual libbio::dispatch_ptr <dispatch_queue_t> consumer_queue() const = 0;
		virtual std::vector <pbwt_sample_type> const &pbwt_samples() const = 0;		// Reduced PBWT samples.
		virtual segmentation_traceback_vector const &reduced_traceback() const = 0;
		virtual permutation_matrix &permutations() = 0;
		virtual bipartite_set_scoring bipartite_set_scoring_method() const = 0;
	
		virtual std::uint32_t sequence_count() const = 0;
		virtual std::uint32_t max_segment_size() const = 0;
	};
}

#endif
