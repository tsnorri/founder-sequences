/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <experimental/iterator>
#include <founder_sequences/segment_text.hh>
#include <libbio/cxxcompat.hh>


namespace founder_sequences {
	
	void segment_text::write_text(
		std::ostream &os,
		segment_text_vector const &slice,
		std::size_t const pos,
		std::size_t const length,
		sequence_vector const &sequences
	) const
	{
		if (is_copied())
			slice[copied_from].write_text(os, slice, pos, length, sequences);
		else
		{
			auto const seq_idx(first_sequence_index());
			auto const &seq(sequences[seq_idx]);
			auto const subspan(seq.subspan(pos, length));
			os.write(reinterpret_cast <char const *>(subspan.data()), subspan.size());
		}
	}
	
	
	// For statistics.
	void segment_text::write(
		std::ostream &os,
		segment_text_vector const &slice,
		std::size_t const pos,
		std::size_t const length,
		sequence_vector const &sequences
	) const
	{
		write_text(os, slice, pos, length, sequences);
		os << '\t';
		
		std::copy(
			sequence_indices.cbegin(),
			sequence_indices.cend(),
			std::experimental::make_ostream_joiner(os, ", ")
		);
		os << '\t';
		
		if (is_copied())
			os << copied_from;
		else
			os << '-';
	}
}
