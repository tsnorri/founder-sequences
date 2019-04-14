/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <algorithm>
#include <founder_sequences/greedy_matcher.hh>
#include <libbio/assert.hh>
#include <libbio/radix_sort.hh>
#include <list>


namespace lb	= libbio;


namespace {
	
	typedef std::uint32_t						seq_index;
	typedef std::pair <seq_index, seq_index>	seq_index_pair;
	typedef std::vector <seq_index>				seq_index_vector;
	typedef std::vector <seq_index_pair>		seq_occurrence_vector;
	typedef std::vector <std::list <seq_index>> index_list_vector;

	
	struct seq_index_pair_access
	{
		seq_index operator()(seq_index_pair const &pair) { return pair.second; }
	};

	
	template <typename t_permutation, typename t_divergence>
	inline void update_string_mappings(
		std::size_t const seq_count,						// Sequence count
		std::size_t const seg_start_pos,					// Segment start position
		t_permutation const &permutation,					// Current permutation
		t_divergence const &divergence,						// Current divergence
		seq_index &distinct_substrings,						// Out: total number of distinct substrings in the segment
		seq_index_vector &seq_index_mapping,				// Out: map {0, 1, …} to the first string (in a_k order) number that begins a run of unique substrings
		seq_index_vector &seq_inverse_mapping,				// Out: inverse of the above
		seq_index_vector &run_lengths						// Out: map {0, 1, …} to each run length
	)
	{
		std::size_t current_run_length(0);	// Length of the current distinct substring run.
		distinct_substrings = 0;
		for (std::size_t i(0); i < seq_count; ++i)
		{
			auto const string_idx(permutation[i]);
			
			// If the divergence value is greater than the segment start position,
			// the current string (in a_k order) starts a run of distinct substrings in this segment.
			if (seg_start_pos < divergence[i])
			{
				run_lengths[distinct_substrings] = current_run_length;	// run_lengths[0] is supposed to be 0.
				seq_index_mapping[distinct_substrings++] = string_idx;	// Map the current value from N to the string number.
				current_run_length = 0;
			}
			
			// Update the inverse mapping.
			libbio_assert(distinct_substrings);
			seq_inverse_mapping[string_idx] = distinct_substrings - 1;
			
			// Increment the run length.
			++current_run_length;
		}
		
		// Store the length of the last run.
		run_lengths[distinct_substrings] = current_run_length;
	}
	
	
	inline void update_seq_occurrences(seq_index const distinct_substrings, seq_index_vector const &run_lengths, seq_occurrence_vector &seq_occurrences)
	{
		// Copy the run lengths to the pairs in seq_occurrences with the sequence indices ({0, 1, …}).
		libbio_assert_eq(0, run_lengths[0]); // run_lengths has 1-based indexing.
		seq_occurrences.clear();
		for (std::size_t i(0); i < distinct_substrings; ++i)
			seq_occurrences.emplace_back(i, run_lengths[1 + i]);
	}
	
	
	void update_copies(
		seq_occurrence_vector const &seq_occurrences,
		std::size_t const max_segment_size,
		std::size_t const seq_count,
		seq_index_vector &cn
	)
	{
		// Fill the copy number vector.
		std::fill(cn.begin(), cn.end(), 1);
		std::size_t const to_fill(max_segment_size - seq_occurrences.size());
		std::size_t rem_size(to_fill);
		
		while (true)
		{
			for (auto const [seq_idx, count] : seq_occurrences)
			{
				std::size_t const copy_count(std::min(rem_size, std::size_t(std::ceil(1.0 * count / seq_count * to_fill))));
				rem_size -= copy_count;
				cn[seq_idx] += copy_count;
				
				if (0 == rem_size)
					return;
			}
		}
	}
	
	
	template <typename t_permutation>
	void update_initial_permutation(
		std::size_t const distinct_substring_count,
		seq_index_vector const &copy_numbers,
		seq_index_vector const &seq_mapping,
		t_permutation &permutation,
		index_list_vector &permutation_slots
	)
	{
		// Copy the (copied) sequence numbers to the permutation and update the slot list.
		libbio_assert_lte(distinct_substring_count, copy_numbers.size());
		std::size_t seq_idx(0);
		std::size_t i(0);
		std::size_t const count(permutation.size());
		for (auto const copy_count : copy_numbers | ranges::view::take_exactly(distinct_substring_count))
		{
			// Update the slot list.
			permutation_slots[seq_idx].clear();
			for (std::size_t j(0); j < copy_count; ++j)
				permutation_slots[seq_idx].emplace_back(j + i);
			
			// Update the permutation.
			std::fill(permutation.begin() + i, permutation.begin() + i + copy_count, seq_mapping[seq_idx]);
			i += copy_count;
			libbio_assert_lte(i, count);
			++seq_idx;
		}
	}
	
	
	inline seq_index find_slot(seq_index const seq_idx, index_list_vector &permutation_slots)
	{
		// Get an available slot from the list.
		// permutation_slots maps {0, 1, 2, …} to lists of permutation indices that have the substring in question.
		auto &slot_list(permutation_slots[seq_idx]);
		libbio_assert(!slot_list.empty());
		auto const slot(slot_list.front());
		slot_list.pop_front();
		return slot;
	}
	
	
	template <typename t_permutation>
	void draw_edge(
		seq_index const lhs_idx,
		seq_index const rhs_idx,
		seq_index_vector const &rhs_seq_mapping,
		index_list_vector &lhs_permutation_slots,
		index_list_vector &rhs_permutation_slots,
		t_permutation &permutation
	)
	{
		// Get a slot.
		auto const slot(find_slot(lhs_idx, lhs_permutation_slots));
		
		// Make the slot available on the right hand side.
		rhs_permutation_slots[rhs_idx].push_back(slot);
		
		// Update the permutation.
		permutation[slot] = rhs_seq_mapping[rhs_idx];
	}
}


namespace founder_sequences
{
	void greedy_matcher::match()
	{
		// Use the greedy algorithm to generate the permutations.
		auto &permutations(m_delegate->permutations()); // Target permutations.
		auto const &pbwt_samples(m_delegate->pbwt_samples());
		
		auto const max_segment_size(m_delegate->max_segment_size());
		assert(max_segment_size);
		auto const matching_bits_needed(lb::bits::highest_bit_set(max_segment_size)); // Space for one additional value.
		auto const matching_max((1 << matching_bits_needed) - 1);
		auto const seq_count(m_delegate->sequence_count());
		libbio_always_assert_lte_msg(2 * matching_bits_needed, 64, "Matching needs to fit in 64 bits");
		
		seq_index lhs_distinct_substrings{};									// Count of distinct substrings in a segment. (Swapped in the main loop below.)
		seq_index rhs_distinct_substrings{};									// Same as above.
		seq_index_vector lhs_seq_mapping(max_segment_size, UINT32_MAX);			// Map {0, 1, 2, …} to the first string (in a_k order) number that begins a run of a unique substring. (Swapped.)
		seq_index_vector rhs_seq_mapping(max_segment_size, UINT32_MAX);			// Same as above.
		seq_index_vector lhs_inverse_seq_mapping(seq_count, UINT32_MAX);		// Inverse of the above but map all string indices in an equivalence class to the same element of N. (Swapped.)
		seq_index_vector rhs_inverse_seq_mapping(seq_count, UINT32_MAX);		// Same as above.
		seq_index_vector lhs_rl(1 + max_segment_size, 0);						// Map {0, 1, 2, …} to the length of the run in question. (Swapped.)
		seq_index_vector rhs_rl(1 + max_segment_size, 0);						// Same as above.
		seq_index_vector lhs_cn(max_segment_size, 0);							// Map {0, 1, 2, …} to the number of copies made in the permutation. (Swapped.)
		seq_index_vector rhs_cn(max_segment_size, 0);							// Same as above.
		seq_index_vector rhs_rc(max_segment_size, 0);							// Copy of the above for bookkeeping (“remaining copies”).
		seq_occurrence_vector seq_occurrences(max_segment_size, {0, 0});		// For sorting {0, 1, 2, …} by the occurrences in the original set.
		seq_occurrence_vector seq_occurrence_buffer(max_segment_size, {0, 0});	// For sorting {0, 1, 2, …} by the occurrences in the original set.
		std::vector <std::uint64_t> index_pairs(seq_count, 0);					// Pairs of matched string indices encoded in std::uint64_t.
		std::vector <std::uint64_t> index_pairs_buffer(seq_count, 0);			// For radix sorting.
		std::size_t seg_start_idx(0);											// Segment start index.
		
		typedef std::list <seq_index_pair> index_pair_list;
		typedef std::map <std::size_t, index_pair_list, std::greater <seq_index>> index_pair_map;
		index_list_vector lhs_permutation_slots(max_segment_size);				// Map {0, 1, 2, …} to permutation slots that have the string in question. (Swapped.)
		index_list_vector rhs_permutation_slots(max_segment_size);				// Same as above.
		index_pair_map index_pairs_by_count;									// Map counts to lists of edges (as pairs of the mapping keys {0, 1, 2, …} above).

		// Fill the string mappings
		{
			auto const &pbwt_sample(pbwt_samples.front());
			auto const &permutation(pbwt_sample.input_permutation());
			auto const &divergence(pbwt_sample.input_divergence());
			
			update_string_mappings(
				seq_count,
				seg_start_idx,
				permutation,
				divergence,
				lhs_distinct_substrings,
				lhs_seq_mapping,
				lhs_inverse_seq_mapping,
				lhs_rl
			);
			
			seg_start_idx = pbwt_sample.sequence_idx();
			update_seq_occurrences(lhs_distinct_substrings, lhs_rl, seq_occurrences);
			lb::radix_sort <true>::sort(seq_occurrences, seq_occurrence_buffer, seq_index_pair_access());
			update_copies(seq_occurrences, max_segment_size, seq_count, lhs_cn);
			update_initial_permutation(lhs_distinct_substrings, lhs_cn, lhs_seq_mapping, permutations.front(), lhs_permutation_slots);
		}

		// Iterate over the PBWT samples. Use a window function in order to get the segment start position of the pair.
		std::size_t target_permutation_idx(1);
		for (auto const &sample : pbwt_samples | ranges::view::drop(1))
		{
			auto const &permutation(sample.input_permutation());
			index_pairs_by_count.clear();
			
			// Fill the string mappings for the right side.
			update_string_mappings(
				seq_count,
				seg_start_idx,
				permutation,
				sample.input_divergence(),
				rhs_distinct_substrings,
				rhs_seq_mapping,
				rhs_inverse_seq_mapping,
				rhs_rl
			);
			
			// Update the distinct subsequence copy numbers.
			update_seq_occurrences(rhs_distinct_substrings, rhs_rl, seq_occurrences);
			lb::radix_sort <true>::sort(seq_occurrences, seq_occurrence_buffer, seq_index_pair_access());
			update_copies(seq_occurrences, max_segment_size, seq_count, rhs_cn);
			
			// Fill the list of edges as encoded index pairs.
			std::size_t current_run_length(0);
			for (std::size_t i(0); i < seq_count; ++i)
			{
				auto const seq_idx(permutation[i]);
				libbio_assert_lt(seq_idx, lhs_inverse_seq_mapping.size());
				libbio_assert_lt(seq_idx, rhs_inverse_seq_mapping.size());
				auto const lhs_idx(lhs_inverse_seq_mapping[seq_idx]);
				auto const rhs_idx(rhs_inverse_seq_mapping[seq_idx]);
				libbio_assert_lte(lhs_idx, matching_max);
				libbio_assert_lte(rhs_idx, matching_max);
				index_pairs[i] = (lhs_idx << matching_bits_needed) | rhs_idx;
			}
			
			// Radix sort the encoded index pairs for counting.
			lb::radix_sort <>::sort(index_pairs, index_pairs_buffer, 2 * matching_bits_needed);
			
			// Map the counts to unique edges. In the worst case, the number of unique counts
			// (index_pair_map keys) is the inverse of 1 + 2 + 3 + … + k + l = k(k + 1) / 2 + l, l ∈ {1, …, k}.
			// Call this w; then w log w = o(m). Hence a balanced binary tree (std::map) may be used.
			{
				std::uint64_t prev_item(index_pairs.front());
				std::size_t current_count(1);				// Count of the current index pair.
				for (std::size_t i(1); i < seq_count; ++i)
				{
					auto const current_item(index_pairs[i]);
					if (prev_item == current_item)
						++current_count;
					else
					{
						// Store the current run of items.
						auto const lhs_idx(prev_item >> matching_bits_needed);
						auto const rhs_idx(prev_item & matching_max);
						libbio_assert_eq(0, (~matching_max) & lhs_idx);
						index_pairs_by_count[current_count].emplace_back(lhs_idx, rhs_idx);
						
						prev_item = current_item;
						current_count = 1;
					}
				}
				{
					// Store the last item.
					auto const lhs_idx(prev_item >> matching_bits_needed);
					auto const rhs_idx(prev_item & matching_max);
					libbio_assert_eq(0, (~matching_max) & lhs_idx);
					index_pairs_by_count[current_count].emplace_back(lhs_idx, rhs_idx);
				}
			}
			
			// Draw the edges.
			{
				// Reset the remaining copy numbers.
				rhs_rc = rhs_cn;
				
				auto &permutation(permutations[target_permutation_idx]);
				bool did_draw_edge(true);
				while (did_draw_edge)
				{
					did_draw_edge = false;
					auto ip_it(index_pairs_by_count.begin());
					auto const ip_end(index_pairs_by_count.end());
					while (ip_it != ip_end)
					{
						// The list below contains the edges that have the current count of occurrences.
						// Try to draw an edge for each.
						auto &list(ip_it->second);
						auto it(list.begin());
						auto const end(list.end());
						while (it != end)
						{
							// Check if there are available substring copies to draw this edge.
							auto const [lhs_idx, rhs_idx] = *it;
							if (lhs_cn[lhs_idx] && rhs_rc[rhs_idx])
							{
								did_draw_edge = true;
								
								// lhs_cn is no longer needed so it may be modified directly.
								--lhs_cn[lhs_idx];
								--rhs_rc[rhs_idx];
								
								draw_edge(lhs_idx, rhs_idx, rhs_seq_mapping, lhs_permutation_slots, rhs_permutation_slots, permutation);
								++it;
							}
							else
							{
								// If not, remove the edge entry.
								it = list.erase(it);
							}
						}
						
						if (!list.empty())
							++ip_it;
						else
						{
							// Remove the list if it was emptied.
							ip_it = index_pairs_by_count.erase(ip_it);
						}
					}
				}
				libbio_assert(index_pairs_by_count.empty());
				
				// Draw the remaining edges.
				{
					std::size_t i(0), j(0);
					while (true)
					{
						// Find available copies on the left side.
						while (0 == lhs_cn[i])
						{
							++i;
							if (i == lhs_distinct_substrings)
								goto end;
						}
						
						// Find available copies on the right side.
						while (0 == rhs_rc[j])
						{
							++j;
							if (j == rhs_distinct_substrings)
								goto end;
						}
						
						draw_edge(i, j, rhs_seq_mapping, lhs_permutation_slots, rhs_permutation_slots, permutation);
						
						--lhs_cn[i];
						--rhs_rc[j];
					}
				end:
					;
				}
			}
			
			using std::swap;
			swap(lhs_distinct_substrings, rhs_distinct_substrings);
			swap(lhs_seq_mapping, rhs_seq_mapping);
			swap(lhs_rl, rhs_rl);
			swap(lhs_cn, rhs_cn);
			swap(lhs_permutation_slots, rhs_permutation_slots);
			seg_start_idx = sample.sequence_idx();

			++target_permutation_idx;
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
