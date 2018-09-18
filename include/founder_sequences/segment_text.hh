/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef FOUNDER_SEQUENCES_SEGMENT_TEXT_HH
#define FOUNDER_SEQUENCES_SEGMENT_TEXT_HH

#include <founder_sequences/founder_sequences.hh>
#include <ostream>
#include <string>
#include <vector>


namespace founder_sequences {
	
	struct segment_text;
	
	typedef std::vector <segment_text>			segment_text_vector;
	typedef std::vector <segment_text_vector>	segment_text_matrix;
	
	
	struct segment_text final
	{
		std::vector <std::size_t>	sequence_indices;
		std::size_t					copied_from{SIZE_MAX};
		
		segment_text() = default;
		
		segment_text(std::size_t copied_from_): copied_from(copied_from_) {}
		
		bool is_copied() const { return SIZE_MAX != copied_from; }
		std::size_t first_sequence_index() const { return sequence_indices.front(); }
		std::uint32_t sequence_count() const { return sequence_indices.size(); }
		std::size_t row_number(std::size_t const row) const { return (SIZE_MAX == copied_from ? row : copied_from); }
		
		// Write the text to the stream.
		void write_text(
			std::ostream &os,
			segment_text_vector const &slice,
			std::size_t const pos,
			std::size_t const length,
			sequence_vector const &sequences
		) const;
		
		// For statistics.
		void write(
			std::ostream &os,
			segment_text_vector const &slice,
			std::size_t const pos,
			std::size_t const length,
			sequence_vector const &sequences
		) const;
	};
}

#endif
