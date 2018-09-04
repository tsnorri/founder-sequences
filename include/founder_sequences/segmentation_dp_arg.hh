/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef FOUNDER_SEQUENCES_SEGMENTATION_DP_ARG_HH
#define FOUNDER_SEQUENCES_SEGMENTATION_DP_ARG_HH

#include <founder_sequences/rmq.hh>
#include <founder_sequences/segment_text.hh>
#include <ostream>
#include <vector>


namespace founder_sequences {
	
	struct segmentation_dp_arg
	{
		std::size_t		lb{0};	// Inclusive.
		std::size_t 	rb{0};	// Exclusive, i.e. [ ) style range.
		std::uint32_t	segment_max_size{UINT32_MAX}; // Maximum of the current segment and the previous one.
		std::uint32_t	segment_size{UINT32_MAX};

		segmentation_dp_arg() = default;
		
		segmentation_dp_arg(
			std::size_t const lb_,
			std::size_t const rb_,
			std::uint32_t const segment_size_
		):
			lb(lb_),
			rb(rb_),
			segment_size(segment_size_)
		{
			assert(lb <= rb);
		}

		segmentation_dp_arg(
			std::size_t const lb_,
			std::size_t const rb_,
			std::uint32_t const segment_max_size_,
			std::uint32_t const segment_size_
		):
			lb(lb_),
			rb(rb_),
			segment_max_size(segment_max_size_),
			segment_size(segment_size_)
		{
			assert(lb <= rb);
		}
		
		bool operator<(segmentation_dp_arg const &rhs) const { return segment_max_size < rhs.segment_max_size; }
		std::size_t text_length() const { return rb - lb; }
	};
	
	std::ostream &operator<<(std::ostream &os, segmentation_dp_arg const &dp_arg);
	
	
	typedef std::vector <segmentation_dp_arg>						segmentation_traceback_vector;
	typedef rmq <segmentation_traceback_vector, std::less <>, 64>	segmentation_traceback_vector_rmq;
	
	
	void output_segments(
		std::ostream &stream,
		segmentation_traceback_vector const &segmentation_traceback,
		segment_text_matrix const &segment_texts
	);
}

#endif
