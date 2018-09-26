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
		substring_copy_number_matrix const &substring_copy_numbers,
		sequence_vector const &sequences
	)
	{
		// Output format:
		// 1. Segment number
		// 2. Left bound (inclusive)
		// 3. Right bound (exclusive)
		// 4. Segment size
		// 5. Subsequence within the segment (string index)
		// 6. Copy number.
		stream << "SEGMENT" "\t" "LB" "\t" "RB" "\t" "SIZE" "\t" "SUBSEQUENCE_NUMBER" "\t" "COPY_NUMBER" "\t" "SUBSEQUENCE" "\n";
		std::size_t segment_idx(0);
		for (auto const &tup : ranges::view::zip(segmentation_traceback, substring_copy_numbers))
		{
			auto const &traceback_arg(std::get <0>(tup));
			auto const &copy_numbers(std::get <1>(tup));
			
			auto const lb(traceback_arg.lb);
			auto const rb(traceback_arg.rb);
			std::size_t prev_copy_number(0);
			
			for (auto const &cn : copy_numbers)
			{
				auto const substring_idx(cn.substring_idx);
				stream << segment_idx << '\t' << lb << '\t' << rb << '\t' << traceback_arg.segment_size << '\t' << substring_idx << '\t' << cn.copy_number - prev_copy_number << '\t';
				prev_copy_number = cn.copy_number;
				
				auto const length(rb - lb);
				auto const &sequence(sequences[substring_idx]);
				auto const &subspan(sequence.subspan(lb, length));
				stream.write(reinterpret_cast <char const *>(subspan.data()), subspan.size());
				
				stream << '\n';
			}
			
			++segment_idx;
		}
		
		stream << std::flush;
	}
	
	
	void output_segments(
		std::ostream &stream,
		segmentation_traceback_vector const &segmentation_traceback,
		segment_text_matrix const &segment_texts,
		sequence_vector const &sequences
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
		for (auto const &tup : ranges::view::zip(segmentation_traceback, segment_texts))
		{
			auto const &traceback_arg(std::get <0>(tup));
			auto const &segment_texts(std::get <1>(tup));
			
			auto const lb(traceback_arg.lb);
			auto const rb(traceback_arg.rb);
			
			for (auto const &seg_text : segment_texts)
			{
				stream << segment_idx << '\t' << lb << '\t' << rb << '\t' << traceback_arg.segment_size << '\t';
				
				seg_text.write_text(
					stream,
					segment_texts,
					lb,
					rb - lb,
					sequences
				);
				stream << '\t';
				
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
