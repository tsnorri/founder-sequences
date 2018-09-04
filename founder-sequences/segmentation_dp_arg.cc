/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <boost/range/combine.hpp>
#include <experimental/iterator>
#include <founder_sequences/segmentation_dp_arg.hh>


namespace founder_sequences {
	
	void output_segments(
		std::ostream &stream,
		segmentation_traceback_vector const &segmentation_traceback,
		segment_text_matrix const &segment_texts
	)
	{
		// Output format (semi-long form):
		// 1. Segment number
		// 2. Left bound (inclusive)
		// 3. Right bound (exclusive)
		// 4. Segment size
		// 5. Subsequence within the segment
		// 6. List of sequence identifiers in which the subsequence occurs separated by commas.
		// 7. From where the subsequence was copied (for bipartite matching).
		stream << "SEGMENT" "\t" "LB" "\t" "RB" "\t" "SIZE" "\t" "SUBSEQUENCE" "\t" "SEQUENCES" "\t" "COPIED_FROM" "\n";
		std::size_t segment_idx(0);
		for (auto const &tup : boost::combine(segmentation_traceback, segment_texts))
		{
			auto const &traceback_arg(tup.get <0>());
			auto const &segment_texts(tup.get <1>());
		
			for (auto const &seg_text : segment_texts)
			{
				stream << segment_idx << '\t' << traceback_arg.lb << '\t' << traceback_arg.rb << '\t' << traceback_arg.segment_size << '\t' << seg_text.text << '\t';
				std::copy(
					seg_text.sequence_indices.cbegin(),
					seg_text.sequence_indices.cend(),
					std::experimental::make_ostream_joiner(stream, ",")
				);
			
				stream
				<< '\t'
				<< (seg_text.is_copied() ? std::to_string(seg_text.copied_from) : "-")
				<< '\n';
			}
			++segment_idx;
		}
		stream << std::flush;
	}
}
