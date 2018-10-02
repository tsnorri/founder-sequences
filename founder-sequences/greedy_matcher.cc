/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <founder_sequences/greedy_matcher.hh>
#include <libbio/radix_sort.hh>


namespace lb	= libbio;


namespace founder_sequences
{
	// Fill index_pairs, to_lhs_substring, to_rhs_string.
	auto greedy_matcher::create_index_pairs(
		matching_tuple const &lhs,
		matching_tuple const &rhs,
		sdsl::int_vector <0> const &rhs_matching,	// From the previous iteration.
		std::vector <detail::substring_index_pair> &index_pairs,
		std::vector <detail::substring_index_pair> &index_pairs_buffer,
		std::vector <std::uint32_t> &lhs_unused_substring_numbers,
		sdsl::int_vector <0> &to_lhs_substring,
		sdsl::int_vector <0> &to_rhs_string
	) const -> std::pair <std::size_t, std::size_t>
	{
		std::size_t lhs_substring_count(0);
		std::size_t rhs_substring_count(0);
		
		auto const &lhs_cn_vector(std::get <0>(lhs));
		auto const &rhs_cn_vector(std::get <0>(rhs));
		auto const &lhs_sample(std::get <1>(lhs));
		auto const &rhs_sample(std::get <1>(rhs));
		assert(lhs_cn_vector.size());
		assert(rhs_cn_vector.size());
		
		// Fill the pairs of matched string indices, count the occurrence of each pair and sort by it.
		// Also update to_lhs_substring s.t. string numbers can be converted to unique substring numbers
		// and to_rhs_string s.t. unique substring numbers can be converted to string numbers.
		assert(0 < rhs_matching.size());
		auto const seq_count(m_delegate->sequence_count());
		auto const max_segment_size(m_delegate->max_segment_size());
		index_pairs.clear();
		index_pairs.resize(seq_count);	// Fill with the default constructor.
		auto const &lhs_input_permutation(lhs_sample.context.input_permutation());
		auto const &rhs_input_permutation(rhs_sample.context.input_permutation());
		std::size_t copy_number(0);
		std::size_t lhs_substring_number(0);
		std::size_t rhs_substring_number(0);
		std::size_t lhs_substring_number_conv(rhs_matching[0]);
		auto lhs_cn_it(lhs_cn_vector.begin());
		auto rhs_cn_it(rhs_cn_vector.begin());
		auto lhs_input_it(lhs_input_permutation.begin());
		auto rhs_input_it(rhs_input_permutation.begin());
		
		while (true)
		{
			assert(lhs_cn_vector.end() != lhs_cn_it);
			assert(rhs_cn_vector.end() != rhs_cn_it);
			
			auto const &lhs_cn(*lhs_cn_it);
			auto const &rhs_cn(*rhs_cn_it);
			auto const limit(std::min(lhs_cn.copy_number, rhs_cn.copy_number));
			
			while (copy_number < limit)
			{
				assert(lhs_input_permutation.end() != lhs_input_it);
				assert(rhs_input_permutation.end() != rhs_input_it);
				
				auto const lhs_string_number(*lhs_input_it++);
				auto const rhs_string_number(*rhs_input_it++);
				
				// Set the values at string positions in index_pairs.
				assert(lhs_string_number < index_pairs.size());
				assert(rhs_string_number < index_pairs.size());
				index_pairs[lhs_string_number].lhs_idx = lhs_substring_number_conv;
				index_pairs[rhs_string_number].rhs_idx = rhs_substring_number;
				
				// Update to_lhs_substring and to to_rhs_string.
				// Note: to_rhs_string gets maximum indices, not minimum. At this time it doesn’t matter, though.
				assert(lhs_string_number < to_lhs_substring.size());
				assert(rhs_substring_number < to_rhs_string.size());
				to_lhs_substring[lhs_string_number] = lhs_substring_number_conv;
				to_rhs_string[rhs_substring_number] = rhs_string_number;
				
				++copy_number;
			}
			
			if (lhs_input_it == lhs_input_permutation.end())
			{
				assert(rhs_input_it == rhs_input_permutation.end());
				break;
			}
			
			if (copy_number == lhs_cn.copy_number)
			{
				++lhs_substring_number;
				++lhs_cn_it;
				lhs_substring_number_conv = rhs_matching[lhs_substring_number];
			}
			
			if (copy_number == rhs_cn.copy_number)
			{
				++rhs_substring_number;
				++rhs_cn_it;
			}
		}
		
		lhs_substring_count = lhs_substring_number + 1;
		rhs_substring_count = rhs_substring_number + 1;
		
		assert(lhs_cn_vector.end() == ++lhs_cn_it);
		assert(rhs_cn_vector.end() == ++rhs_cn_it);
		assert(lhs_substring_count == lhs_cn_vector.size());
		assert(rhs_substring_count == rhs_cn_vector.size());
		
		assert(rhs_substring_count < to_rhs_string.size());
		std::fill(to_rhs_string.begin() + rhs_substring_count, to_rhs_string.end(), m_permutation_max);
		
		// Store the matching indices that represent the dash sequence.
		lhs_unused_substring_numbers.clear();
		std::copy(
			rhs_matching.begin() + lhs_substring_count,
			rhs_matching.begin() + max_segment_size,
			std::back_inserter(lhs_unused_substring_numbers)
		);
		
		// Set the counts.
		for (auto &pair : index_pairs)
			pair.count = 1;
		
		// Sort by the indices. The sorting algorithm is stable.
		// This might be slow, though, since each value is stored in a separate variable and the highest bit needs to be determined separately for each of them.
		lb::radix_sort <>::sort(index_pairs, index_pairs_buffer, [](detail::substring_index_pair const &index_pair){ return index_pair.rhs_idx; });
		lb::radix_sort <>::sort(index_pairs, index_pairs_buffer, [](detail::substring_index_pair const &index_pair){ return index_pair.lhs_idx; });
		
		// Find equivalent pairs, store counts.
		lb::unique_count(index_pairs.begin(), index_pairs.end(), index_pairs_buffer);
		
		// Swap s.t. index_pairs contains the pairs.
		{
			using std::swap;
			swap(index_pairs, index_pairs_buffer);
		}
		
		// index_pairs should now be sorted by lhs_idx, then rhs_idx.
		assert(std::is_sorted(index_pairs.cbegin(), index_pairs.cend(), [](detail::substring_index_pair const &lhs, detail::substring_index_pair const &rhs) -> bool {
			if (lhs.lhs_idx < rhs.lhs_idx)
				return true;
			else if (rhs.lhs_idx < lhs.lhs_idx)
				return false;
			else if (lhs.rhs_idx < rhs.rhs_idx)
				return true;
			
			return false;
		}));
		
		// Reverse sort by count.
		lb::radix_sort <true>::sort(index_pairs, index_pairs_buffer, [](detail::substring_index_pair const &index_pair){ return index_pair.count; });
		
		return std::make_pair(lhs_substring_count, rhs_substring_count);
	}
	
	
	// Fill lhs_matching and rhs_matching.
	void greedy_matcher::create_matching(
		std::vector <detail::substring_index_pair> const &index_pairs,
		std::uint64_t const matching_max,
		sdsl::int_vector <0> &lhs_matching,
		sdsl::int_vector <0> &rhs_matching
	) const
	{
		auto const max_segment_size(m_delegate->max_segment_size());
		
		// Create the matching by iterating the index pairs in descending occurrence order.
		// FIXME: are both lhs_matching and rhs_matching actually needed?
		std::fill(lhs_matching.begin(), lhs_matching.end(), matching_max);
		std::fill(rhs_matching.begin(), rhs_matching.end(), matching_max);
		for (auto const &pair : index_pairs)
		{
			assert(pair.lhs_idx < lhs_matching.size());
			assert(pair.rhs_idx < rhs_matching.size());
			
			if (lhs_matching[pair.lhs_idx] == matching_max && rhs_matching[pair.rhs_idx] == matching_max)
			{
				lhs_matching[pair.lhs_idx] = pair.rhs_idx;
				rhs_matching[pair.rhs_idx] = pair.lhs_idx;
			}
		}
		
		// Draw the remaining edges.
		std::size_t lhs_idx(0);
		std::size_t rhs_idx(0);
		while (true)
		{
			while (lhs_idx < max_segment_size && lhs_matching[lhs_idx] != matching_max)
				++lhs_idx;
			
			while (rhs_idx < max_segment_size && rhs_matching[rhs_idx] != matching_max)
				++rhs_idx;
			
			if (lhs_idx == max_segment_size || rhs_idx == max_segment_size)
			{
				assert(lhs_idx == rhs_idx);
				break;
			}
			
			lhs_matching[lhs_idx] = rhs_idx;
			rhs_matching[rhs_idx] = lhs_idx;
		}
	}

	
	void greedy_matcher::match()
	{
		// Use the greedy algorithm to generate the permutations.
		auto &permutations(m_delegate->permutations());
		
		auto const max_segment_size(m_delegate->max_segment_size());
		assert(max_segment_size);
		auto const matching_bits_needed(lb::bits::highest_bit_set(max_segment_size)); // Space for one additional value.
		auto const matching_max((1 << matching_bits_needed) - 1);
		auto const seq_count(m_delegate->sequence_count());

		std::vector <detail::substring_index_pair> index_pairs;			// Pairs of matched string indices.
		std::vector <detail::substring_index_pair> index_pairs_buffer;	// For radix sorting and eliminating duplicates.
		std::vector <std::uint32_t> lhs_unused_substring_numbers;
		sdsl::int_vector <0> lhs_matching(max_segment_size, matching_max, matching_bits_needed);
		sdsl::int_vector <0> rhs_matching(max_segment_size, 0, matching_bits_needed);
		sdsl::int_vector <0> to_lhs_substring(seq_count, 0, matching_bits_needed);		// For converting the string numbers to unique substring numbers.
		sdsl::int_vector <0> to_rhs_string(1 + matching_max, 0, m_permutation_bits_needed);

		// Set the initial state with the substring numbers and the minimum unique substring indices.
		// rhs_matching is used first.
		std::iota(rhs_matching.begin(), rhs_matching.end(), 0);	// Start with identity.
		{
			auto const &cn_vec(m_substrings_to_output->front());
			auto &permutation(permutations.front());
			
			std::size_t i(0);
			auto const cn_vec_size(cn_vec.size());
			while (i < cn_vec_size)
			{
				auto const &cn(cn_vec[i]);
				auto const substring_idx(cn.substring_idx);
				permutation[i] = substring_idx;
				++i;
			}
			
			// Assign numbers for gap sequences.
			auto const permutation_size(permutation.size());
			while (i < permutation_size)
			{
				permutation[i] = m_permutation_max;
				++i;
			}
		}
		
		for (auto const &tup : ranges::view::zip(*m_substrings_to_output, m_delegate->pbwt_samples(), permutations) | ranges::view::sliding(2))
		{
			auto const &lhs(tup[0]);
			auto const &rhs(tup[1]);
			auto const &lhs_output_permutation(std::get <2>(lhs));
			auto &rhs_output_permutation(std::get <2>(rhs));
			
			// Create pairs of string indices of the unique substrings and sort them.
			auto const res(create_index_pairs(
				lhs,
				rhs,
				rhs_matching,
				index_pairs,
				index_pairs_buffer,
				lhs_unused_substring_numbers,
				to_lhs_substring,
				to_rhs_string
			));
			std::size_t lhs_substring_count(res.first);
			std::size_t rhs_substring_count(res.second);
			libbio_assert_eq(max_segment_size - lhs_substring_count, lhs_unused_substring_numbers.size());
			create_matching(index_pairs, matching_max, lhs_matching, rhs_matching);
			
			{
				// Create the next permutation.
				std::size_t idx(0);
				auto lhs_unused_it(lhs_unused_substring_numbers.cbegin());
				for (auto const lhs_string_idx : lhs_output_permutation)
				{
					// Check if there is a gap sequence on the left.
					// Join to the gap sequences in an arbitrary order.
					std::size_t lhs_substring_idx(0);
					if (lhs_string_idx == m_permutation_max)
						lhs_substring_idx = *lhs_unused_it++;
					else
						lhs_substring_idx = to_lhs_substring[lhs_string_idx];
					
					auto const matched_idx(lhs_matching[lhs_substring_idx]);
					auto const rhs_string_idx(to_rhs_string[matched_idx]);
					rhs_output_permutation[idx++] = rhs_string_idx;
				}
				
				assert(lhs_unused_it == lhs_unused_substring_numbers.cend());
			}
		}
		
		m_delegate->matcher_did_finish(*this);
	}
	
	
	void greedy_matcher::output_segments(std::ostream &stream, sequence_vector const &sequences)
	{
		::founder_sequences::output_segments(
			stream,
			m_delegate->reduced_traceback(),
			*m_substrings_to_output,
			sequences
		);
	}
}
