/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <experimental/iterator>
#include <founder_sequences/segment_text.hh>


namespace founder_sequences {
	
	std::ostream &operator<<(std::ostream &os, segment_text const &seg_text)
	{
		os << seg_text.text << '\t';
		std::copy(
			seg_text.sequence_indices.cbegin(),
			seg_text.sequence_indices.cend(),
			std::experimental::make_ostream_joiner(os, ", ")
		);
		os << '\t';
		
		if (seg_text.is_copied())
			os << seg_text.copied_from;
		else
			os << '-';
		
		return os;
	}
}
