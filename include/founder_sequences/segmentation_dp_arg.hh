/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef FOUNDER_SEQUENCES_SEGMENTATION_DP_ARG_HH
#define FOUNDER_SEQUENCES_SEGMENTATION_DP_ARG_HH

#include <founder_sequences/rmq.hh>
#include <founder_sequences/segment_text.hh>
#include <founder_sequences/substring_copy_number.hh>
#include <ostream>
#include <vector>


namespace founder_sequences {
	
	struct segmentation_dp_arg final
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
	
	inline std::ostream &operator<<(std::ostream &os, segmentation_dp_arg const &dp_arg)
	{
		os << "([" << dp_arg.lb << ", " << dp_arg.rb << ") s: " << dp_arg.segment_size << " ms: " << dp_arg.segment_max_size << ")";
		return os;
	}
	
	
	typedef std::vector <segmentation_dp_arg>						segmentation_traceback_vector;
	typedef rmq <segmentation_traceback_vector, std::less <>, 64>	segmentation_traceback_vector_rmq;
	
	void output_segments(
		std::ostream &stream,
		segmentation_traceback_vector const &segmentation_traceback,
		substring_copy_number_matrix const &substring_copy_numbers,
		sequence_vector const &sequences
	);
	
	void output_segments(
		std::ostream &stream,
		segmentation_traceback_vector const &segmentation_traceback,
		segment_text_matrix const &segment_texts,
		sequence_vector const &sequences
	);
}


namespace boost { namespace serialization {
	
	template <typename t_archive>
	void serialize(t_archive &ar, founder_sequences::segmentation_dp_arg &arg, unsigned int const version)
	{
		ar & arg.lb;
		ar & arg.rb;
		ar & arg.segment_max_size;
		ar & arg.segment_size;
	}
}}

#endif
