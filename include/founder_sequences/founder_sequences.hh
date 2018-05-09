/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef FOUNDER_SEQUENCES_FOUNDER_SEQUENCES_HH
#define FOUNDER_SEQUENCES_FOUNDER_SEQUENCES_HH

#include <string>
#include <vector>


namespace founder_sequences {
	
	struct segment_text
	{
		std::string					text;
		std::vector <std::size_t>	sequence_indices;
		std::size_t					copied_from{SIZE_MAX};
		
		bool is_copied() const { return SIZE_MAX != copied_from; }
		std::size_t first_sequence_index() const { return sequence_indices.front(); }
		std::size_t sequence_count() const { return sequence_indices.size(); }
		std::size_t row_number(std::size_t const row) const { return (SIZE_MAX == copied_from ? row : copied_from); }
	};
	
	std::ostream &operator<<(std::ostream &os, segment_text const &seg_text);
	
	typedef std::vector <segment_text>			segment_text_vector;
	typedef std::vector <segment_text_vector>	segment_text_matrix;
	
	
	typedef std::vector <std::uint32_t>			matching_vector;
	
	
	enum class segment_joining : uint8_t {
		MATCHING	= 0,
		RANDOM
	};
	
	enum class input_format : uint8_t {
		FASTA		= 0,
		LIST_FILE
	};
}

#endif
