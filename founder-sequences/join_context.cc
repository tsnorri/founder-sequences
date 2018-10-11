/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <founder_sequences/join_context.hh>
#include <libbio/algorithm.hh>
#include <libbio/bits.hh>

namespace lb = libbio;


namespace founder_sequences {

	void join_context::make_cumulative_sum(substring_copy_number_vector &vec) const
	{
		std::size_t cumulative_sum(0);
		for (auto &cn : vec)
		{
			cumulative_sum += cn.copy_number;
			cn.copy_number = cumulative_sum;
		}
	}
	
	
	void join_context::output_gaps(
		std::ostream &ostream,
		std::size_t const text_pos,
		std::size_t const text_length
	) const
	{
		std::fill_n(std::ostreambuf_iterator <char>(ostream), text_length, '-');
	}
	
	
	void join_context::output_substring(
		std::ostream &ostream,
		std::size_t const sequence_idx,
		std::size_t const text_pos,
		std::size_t const text_length,
		sequence_vector const &sequences
	) const
	{
		auto const sequence(sequences[sequence_idx]);
		auto const subspan(sequence.subspan(text_pos, text_length));
		auto const bytes(std::as_bytes(subspan));
		ostream.write(reinterpret_cast <char const *>(bytes.data()), bytes.size());
	}
	
	
	void join_context::join_segments_and_output(segment_joining const seg_joining)
	{
		assert(dispatch_get_current_queue() == dispatch_get_main_queue());
		
		// Count the instances of each substring.
		assert(m_segmentation_container.reduced_traceback.size() == m_segmentation_container.reduced_pbwt_samples.size());
		m_substring_copy_numbers.clear();
		m_substring_copy_numbers.resize(m_segmentation_container.reduced_traceback.size());
		
		lb::parallel_for_each(
			ranges::view::zip(m_segmentation_container.reduced_traceback, m_segmentation_container.reduced_pbwt_samples, m_substring_copy_numbers),
			[this, seg_joining](auto const &tup){
				auto const &dp_arg(std::get <0>(tup));
				auto const &sample(std::get <1>(tup));
				auto &substring_cn(std::get <2>(tup));
				
				// Count the instances w.r.t. dp_arg’s left bound and sort in decreasing order.
				// Then, in case of non-greedy matching, fill the segment up to the maximum
				// segment size by copying substrings in proportion to their occurrence.
				auto const substring_count(sample.unique_substring_count_idxs_lhs(dp_arg.lb, substring_cn));
				assert(substring_count);
				assert(substring_cn.size());
				
				// Numbering needed for PBWT order matching.
				{
					std::uint32_t i(0);
					for (auto &cn : substring_cn)
						cn.string_idx = i++;
				}
				
				if (! (segment_joining::GREEDY == seg_joining || segment_joining::BIPARTITE_MATCHING == seg_joining))
				{
					// Sort by count.
					std::sort(substring_cn.begin(), substring_cn.end());
					
					// Assign new copy numbers in proportion.
					auto const empty_slots(m_segmentation_container.max_segment_size - substring_count);
					std::size_t remaining_slots(empty_slots);
					for (auto &cn : substring_cn | ranges::view::reverse)
					{
						auto const addition(lb::min_ct(remaining_slots, std::ceil(1.0 * cn.copy_number / substring_count * empty_slots)));
						cn.copy_number = 1 + addition;
						remaining_slots -= addition;
					}
					
					// If there are still slots left, add to copy numbers.
					while (remaining_slots)
					{
						for (auto &cn : substring_cn | ranges::view::reverse)
						{
							++cn.copy_number;
							--remaining_slots;
							if (0 == remaining_slots)
								goto loop_end;
						}
					}
					
				loop_end:
					// Check that the sum of copy numbers matches m_max_segment_size.
					assert(m_segmentation_container.max_segment_size == ranges::accumulate(substring_cn | ranges::view::transform([](auto const &cn) -> std::size_t { return cn.copy_number; }), 0));

					// For PBWT order output sort in the original order.
					if (segment_joining::PBWT_ORDER == seg_joining)
					{
						std::sort(substring_cn.begin(), substring_cn.end(), [](substring_copy_number const &lhs, substring_copy_number const &rhs){
							return lhs.string_idx < rhs.string_idx;
						});
					}
				}
				
				make_cumulative_sum(substring_cn);
			}
		);
		
		libbio_assert_eq(0, std::count_if(
			m_substring_copy_numbers.cbegin(),
			m_substring_copy_numbers.cend(),
			[](auto const &vec) -> bool { return 0 == vec.size(); }
		));
		
		// *this may be invalid after calling context_did_output_founders().
		switch (seg_joining)
		{
			case segment_joining::GREEDY:
				join_greedy();
				break;
				
			case segment_joining::BIPARTITE_MATCHING:
				join_with_bipartite_matching();
				break;
				
			case segment_joining::RANDOM:
				join_random_order_and_output();
				break;
				
			case segment_joining::PBWT_ORDER:
				m_delegate->context_will_output_founders(*this);
				join_pbwt_order_and_output();
				m_delegate->context_did_output_founders(*this);
				break;
				
			default:
				libbio_fail("Unexpected segment joining method.");
				break;
		}
	}
	
	
	void join_context::output_segments(segment_joining const seg_joining) const
	{
		auto &stream(m_delegate->segments_output_stream());
		auto const &sequences(m_delegate->sequences());
		switch (seg_joining)
		{
			case segment_joining::GREEDY:
			case segment_joining::BIPARTITE_MATCHING:
				m_matcher->output_segments(stream, sequences);
				break;
				
			case segment_joining::RANDOM:
			case segment_joining::PBWT_ORDER:
				::founder_sequences::output_segments(
					stream,
					m_segmentation_container.reduced_traceback,
					m_substring_copy_numbers,
					sequences
				);
					break;
				
			default:
				libbio_fail("Unexpected segment joining method.");
				break;
		}
	}
	
	
	void join_context::init_permutations()
	{
		// Since we output row-wise, all of the permutation vectors need to be pre-calculated.
		// Try to save some space by using SDSL’s integer vectors.
		// seq_count is not a string index but the following will allow storing it.
		auto const seq_count(sequence_count());
		m_permutation_bits_needed = lb::bits::highest_bit_set(seq_count);
		m_permutation_max = (1 << m_permutation_bits_needed) - 1;
		assert(seq_count <= m_permutation_max);
		
		m_permutations.clear();
		m_permutations.resize(m_segmentation_container.reduced_traceback.size());
		for (auto &permutation : m_permutations)
		{
			sdsl::int_vector <0> temp_permutation(m_segmentation_container.max_segment_size, 0, m_permutation_bits_needed);
			permutation = std::move(temp_permutation);
		}
	}
	
	
	void join_context::join_greedy()
	{
		init_permutations();
		m_matcher.reset(new greedy_matcher(*this, m_substring_copy_numbers, m_permutation_max, m_permutation_bits_needed));
		m_matcher->match();
	}
	
	
	void join_context::join_with_bipartite_matching()
	{
		init_permutations();
		m_matcher.reset(new bipartite_matcher(*this, m_substring_copy_numbers));
		m_matcher->match();
	}
	
	
	void join_context::matcher_did_finish(bipartite_matcher &matcher)
	{
		assert(dispatch_get_current_queue() == dispatch_get_main_queue());
		
		m_delegate->context_will_output_founders(*this);
		output_in_permutation_order();
		m_delegate->context_did_output_founders(*this);
	}
	
	
	void join_context::matcher_did_finish(greedy_matcher &matcher)
	{
		assert(dispatch_get_current_queue() == dispatch_get_main_queue());
		
		m_delegate->context_will_output_founders(*this);
		output_in_permutation_order();
		m_delegate->context_did_output_founders(*this);
	}
	
	
	void join_context::join_random_order_and_output()
	{
		assert(dispatch_get_current_queue() == dispatch_get_main_queue());
		
		init_permutations();
	
		// Instantiate an std::mt19937 with the given seed and generate a permutation of the indices for each segment.
		std::mt19937 urbg(m_random_seed);
	
		// Fill the permutations and shuffle.
		for (auto const &tup : ranges::view::zip(m_substring_copy_numbers, m_permutations))
		{
			auto const &cn_vector(std::get <0>(tup));
			auto &permutation(std::get <1>(tup));
			auto it(permutation.begin());
			for (auto const &cn : cn_vector)
			{
				auto end(it + cn.copy_number);
				std::fill(it, end, cn.substring_idx);
				it = end;
			}
			std::shuffle(permutation.begin(), permutation.end(), urbg);
		}
		
		m_delegate->context_will_output_founders(*this);
		output_in_permutation_order();
		m_delegate->context_did_output_founders(*this);
	}
	
	
	void join_context::join_pbwt_order_and_output() const
	{
		assert(dispatch_get_current_queue() == dispatch_get_main_queue());
		
		auto const &sequences(m_delegate->sequences());
		auto &stream(m_delegate->sequence_output_stream());
		
		// Create iterators to each vector of runs.
		std::vector <substring_copy_number_vector::const_iterator> cn_iterators(m_substring_copy_numbers.size());
		for (auto const &tup : ranges::view::zip(m_substring_copy_numbers, cn_iterators))
			std::get <1>(tup) = std::get <0>(tup).cbegin();
		
		// Output.
		for (std::size_t row(0); row < m_segmentation_container.max_segment_size; ++row)
		{
			for (auto const &tup : ranges::view::zip(cn_iterators, m_segmentation_container.reduced_traceback))
			{
				auto &it(std::get <0>(tup));
				auto const &dp_arg(std::get <1>(tup));
				auto const ending_index(it->copy_number); // Cumulative sum.
				
				if (row == ending_index)
					++it;
				
				assert(row < it->copy_number);
				
				auto const substring_idx(it->substring_idx);
				output_substring(stream, substring_idx, dp_arg.lb, dp_arg.text_length(), sequences);
			}
			stream << '\n';
		}
		
		stream << std::flush;
		
#ifndef NDEBUG
		for (auto const &tup : ranges::view::zip(m_substring_copy_numbers, cn_iterators))
			assert(++std::get <1>(tup) == std::get <0>(tup).cend());
#endif
	}
	
	
	void join_context::output_in_permutation_order() const
	{
		auto const &sequences(m_delegate->sequences());
		auto &stream(m_delegate->sequence_output_stream());
		
		assert(m_permutations.size() == m_segmentation_container.reduced_traceback.size());
		
		// Output.
		for (std::size_t row(0); row < m_segmentation_container.max_segment_size; ++row)
		{
			for (auto const &tup : ranges::view::zip(m_permutations, m_segmentation_container.reduced_traceback))
			{
				auto &permutation(std::get <0>(tup));
				auto const &dp_arg(std::get <1>(tup));
				auto const substring_idx(permutation[row]);
				if (m_permutation_max == substring_idx)
					output_gaps(stream, dp_arg.lb, dp_arg.text_length());
				else
					output_substring(stream, substring_idx, dp_arg.lb, dp_arg.text_length(), sequences);
			}
			stream << '\n';
		}
		
		stream << std::flush;
	}
}
