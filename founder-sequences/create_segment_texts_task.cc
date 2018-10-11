/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <founder_sequences/create_segment_texts_task.hh>
#include <libbio/algorithm.hh>


namespace lb	= libbio;


namespace founder_sequences {
	
	void create_segment_texts_task::execute()
	{
		m_segment_texts->resize(m_max_segment_size);
		auto const permutation(m_pbwt_sample->input_permutation());
		std::size_t seg_idx(0);
		std::size_t string_idx(0);
		
		// Create a segment_text for each distinct substring.
		for (auto const &cn : *m_substring_copy_numbers)
		{
			segment_text seg;
			seg.sequence_indices.resize(cn.copy_number - string_idx, 0);
			std::size_t i(0);
			while (string_idx < cn.copy_number)
				seg.sequence_indices[i++] = permutation[string_idx++];
			
			// FIXME: use (MSB?) radix sort?
			std::sort(seg.sequence_indices.begin(), seg.sequence_indices.end());
			(*m_segment_texts)[seg_idx++] = std::move(seg);
		}
		
		// If there are remaining slots, fill them.
		if (seg_idx < m_max_segment_size)
		{
			auto const empty_slots(m_max_segment_size - seg_idx);
			std::size_t remaining_slots(empty_slots);
			
			// Sort the items by source text count.
			std::sort(m_segment_texts->begin(), m_segment_texts->begin() + seg_idx, [](segment_text const &lhs, segment_text const &rhs) -> bool {
				return lhs.sequence_count() > rhs.sequence_count();
			});
			
			// Copy the segment_texts if needed.
			auto const limit(seg_idx);
			auto it(m_segment_texts->begin() + seg_idx);
			for (std::size_t i(0); i < limit; ++i)
			{
				auto const &seg((*m_segment_texts)[i]);
				auto const copy_number(lb::min_ct(remaining_slots, std::ceil(1.0 * seg.sequence_count() / m_seq_count * remaining_slots)));
				
				segment_text copied_seg(i);
				it = std::fill_n(it, copy_number, copied_seg);
				
				remaining_slots -= copy_number;
				if (0 == remaining_slots)
					break;

				assert(it != m_segment_texts->end());
			}

			// If there are remaining slots, fill them.
			while (remaining_slots)
			{
				for (std::size_t i(0); i < limit; ++i)
				{
					*it++ = segment_text(i);
					--remaining_slots;
					
					if (0 == remaining_slots)
						goto loop_end;
					
					assert(m_segment_texts->end() != it);
				}
			}
		loop_end:
			assert(it == m_segment_texts->end());
		}
	}
}
