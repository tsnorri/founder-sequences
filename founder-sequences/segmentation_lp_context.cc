/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <founder_sequences/segmentation_lp_context.hh>
#include <libbio/algorithm.hh>
#include <libbio/bits.hh>
#include <libbio/radix_sort.hh>

namespace lb = libbio;


namespace founder_sequences {
	
	void segmentation_lp_context::generate_traceback(std::size_t const lb, std::size_t const rb)
	{
		// Calculate the first L - 1 columns, which gives the required result for calculating M(L).
		// idx is 0-based, m_segment_length is 1-based.
		dispatch_async(*m_producer_queue, ^{
			m_pbwt_ctx.prepare();
		
			auto const seq_length(m_pbwt_ctx.sequence_size());
			auto const seq_count(m_pbwt_ctx.size());
			auto const segment_length(m_delegate->segment_length());
			auto const dp_size(seq_length - segment_length + 1);
		
			{
				// Values shifted to the left by m_segment_length (L) since the first L columns have the same value anyway.
				segmentation_traceback_vector temp(dp_size);
				segmentation_traceback_vector_rmq temp_rmq(temp);
			
				using std::swap;
				swap(m_segmentation_traceback_dp, temp);
				swap(m_segmentation_traceback_dp_rmq, temp_rmq);
				m_segmentation_traceback_dp_rmq.set_values(m_segmentation_traceback_dp);
			}
			
			m_pbwt_ctx.process(
				segment_length - 1,
				[this](std::uint32_t const idx, auto const &){ output_segmentation_status_mq(1 + idx); },
				[this, lb, rb](){
					dispatch_async(*m_producer_queue, ^{
						generate_traceback_part_2(lb, rb);
					});
				}
			);
		});
	}
	
	
	void segmentation_lp_context::generate_traceback_part_2(std::size_t const lb, std::size_t const rb)
	{
		// Calculate the columns L to 2L - 1 (1-based), which gives the size of one segment for each column since the minimum
		// length of two segments is 2L.
		dispatch_async(*m_producer_queue, ^{
			auto const seq_count(m_pbwt_ctx.size());
			auto const segment_length(m_delegate->segment_length());
			auto const limit(std::min(2 * segment_length, rb - segment_length) - 1);
			
			m_pbwt_ctx.process(
				limit,
				[this, lb, seq_count, segment_length](std::uint32_t const idx, auto const &counts){
					// Counts is just a column object, copy.
					dispatch_async(*m_consumer_queue, ^{
						// Calculate the segment size by finding the range of the relevant key
						// (which is “0” in this case b.c. the column count is less than 2L,
						// so the range is at most [counts.begin(), counts.begin() + 1))
						// and subtracting the count from the sequence count.
						std::size_t segment_size_diff(0);
						auto const begin(counts.cbegin_pairs());
						if (counts.cend_pairs() != begin && 0 == begin->first)
							segment_size_diff = begin->second;
					
						auto const tb_idx(idx + 1 - segment_length);
						auto const segment_size(seq_count - segment_size_diff);
						segmentation_dp_arg const current_arg(lb, 1 + idx, segment_size, segment_size);
						m_segmentation_traceback_dp[tb_idx] = current_arg;
						m_segmentation_traceback_dp_rmq.update(tb_idx);
					
						output_segmentation_status_mq(1 + idx);
					});
				},
				[this, lb, rb](){
					dispatch_async(*m_producer_queue, ^{
						generate_traceback_part_3(lb, rb);
					});
				}
			);
		});
	}
	
	
	void segmentation_lp_context::generate_traceback_part_3(std::size_t const lb, std::size_t const rb)
	{
		// Calculate the columns (L or) 2L to n (1-based). This results in multiple segments the minimum length of each is L.
		// Both rb an m_segment_length are 1-based.
		dispatch_async(*m_producer_queue, ^{
			auto const seq_count(m_pbwt_ctx.size());
			auto const segment_length(m_delegate->segment_length());
			auto const limit(rb - segment_length);
			
			m_pbwt_ctx.process(
				limit,
				[this, lb, seq_count, segment_length](std::uint32_t const idx, auto const &counts){
					dispatch_async(*m_consumer_queue, ^{
						// Use the texts up to this point as the initial value.
						segmentation_dp_arg min_arg(lb, 1 + idx, seq_count, seq_count);
						calculate_segmentation_lp_dp_arg(
							counts,
							m_segmentation_traceback_dp,
							m_segmentation_traceback_dp_rmq,
							seq_count,
							segment_length,
							lb,
							idx,
							min_arg
						);
					
						auto const tb_idx(idx + 1 - segment_length);
						m_segmentation_traceback_dp[tb_idx] = min_arg;
						m_segmentation_traceback_dp_rmq.update(tb_idx);
					
						output_segmentation_status_mq(1 + idx);
					});
				},
				[this, lb, rb](){
					dispatch_async(*m_producer_queue, ^{
						generate_traceback_part_4(lb, rb);
					});
				}
			);
		});
	}
	
	
	void segmentation_lp_context::generate_traceback_part_4(std::size_t const lb, std::size_t const rb)
	{
		// Calculate the size of the final segment.
		dispatch_async(*m_producer_queue, ^{
			auto const seq_count(m_pbwt_ctx.size());
			auto const segment_length(m_delegate->segment_length());
			
			m_pbwt_ctx.process(
				rb,
				[this](std::uint32_t const idx, auto const &){
					output_segmentation_status_mq(1 + idx);
				},
				[this, lb, rb, seq_count, segment_length](){
					dispatch_async(dispatch_get_main_queue(), ^{
						std::cerr << std::endl;
					});
					
					dispatch_async(*m_consumer_queue, ^{
						auto counts(m_pbwt_ctx.last_divergence_value_counts());
						auto const idx(m_pbwt_ctx.sequence_idx());
						segmentation_dp_arg min_arg(lb, idx, seq_count, seq_count);
						calculate_segmentation_lp_dp_arg(
							counts,
							m_segmentation_traceback_dp,
							m_segmentation_traceback_dp_rmq,
							seq_count,
							segment_length,
							lb,
							idx - 1,
							min_arg
						);
					
						// Position of the last traceback argument.
						auto const tb_idx(m_segmentation_traceback_dp.size() - 1);
						assert(rb - segment_length == tb_idx);
					
						m_segmentation_traceback_dp[tb_idx] = min_arg;
						
						follow_traceback();
					});
				}
			);
		});
	}
	
	
	void segmentation_lp_context::follow_traceback()
	{
		// Follow the traceback.
		dispatch_async(dispatch_get_main_queue(), ^{
			std::cerr << "Following the traceback…" << std::endl;
		});
		
		assert(m_segmentation_traceback_dp.size());
		m_segmentation_traceback_res.clear();
		std::size_t arg_idx(m_segmentation_traceback_dp.size() - 1);
		auto const segment_length(m_delegate->segment_length());
		while (true)
		{
			auto const &current_arg(m_segmentation_traceback_dp[arg_idx]);
			// Fill in reverse order.
			m_segmentation_traceback_res.push_back(current_arg);
		
			auto const next_pos(current_arg.lb);
			if (0 == next_pos) // FIXME: generalize by substituting 0 with lb.
				break;
		
			assert(segment_length <= next_pos);
			arg_idx = next_pos - segment_length;
		}
	
		// Reverse the filled traceback.
		std::reverse(m_segmentation_traceback_res.begin(), m_segmentation_traceback_res.end());
		
		// Store the maximum size.
		m_max_segment_size = m_segmentation_traceback_res.back().segment_max_size;
		
		dispatch_async(dispatch_get_main_queue(), ^{
			m_delegate->context_did_finish_traceback(*this);
		});
	}
	
	
	// Given a collection of PBWT samples, update them s.t. their right bound matches that of
	// the traceback arguments.
	void segmentation_lp_context::update_samples_to_traceback_positions()
	{
		m_update_samples_group.reset(dispatch_group_create());
		
		auto &pbwt_samples(m_pbwt_ctx.samples());
		auto const sample_count(pbwt_samples.size());
		auto traceback_it(m_segmentation_traceback_res.cbegin());
		auto const traceback_end(m_segmentation_traceback_res.cend());
		
		std::size_t lb(0);
		std::size_t i(1);
		std::size_t last_moved_sample(SIZE_MAX);
		text_position_vector right_bounds;
		while (i < sample_count)
		{
			if (traceback_it == traceback_end)
				break;
			
			lb = traceback_it->lb;
			
			// Find the traceback arguments that are located between the previous sample (i - 1) and the current one.
			auto const &next_sample(pbwt_samples[i]);
			while (traceback_it < traceback_end && traceback_it->rb < next_sample.rb)
			{
				right_bounds.emplace_back(traceback_it->rb);
				++traceback_it;
			}
			
			// If any were found, start a task.
			assert(ranges::is_sorted(right_bounds));
			if (right_bounds.size())
			{
				last_moved_sample = i - 1;
				auto &prev_sample(pbwt_samples[last_moved_sample]);
				
				// Take the accumulated right_bounds.
				text_position_vector current_right_bounds;
				using std::swap;
				swap(right_bounds, current_right_bounds);
				
				// Start the task.
				start_update_sample_task(lb, std::move(prev_sample), std::move(current_right_bounds));
			}
			
			++i;
		}
		
		// Handle the remaining traceback arguments.
		{
			std::transform(traceback_it, traceback_end, std::back_inserter(right_bounds), [](auto const &tb) {
				return tb.rb;
			});

			assert(ranges::is_sorted(right_bounds));
			if (right_bounds.size())
			{
				assert(i - 1 != last_moved_sample);
				auto &last_sample(pbwt_samples[i - 1]);
				assert(last_sample.rb < right_bounds.front());
				start_update_sample_task(lb, std::move(last_sample), std::move(right_bounds));
			}
		}
		
		// Remove the remaining samples.
		pbwt_samples.clear();
		
		dispatch_group_notify(*m_update_samples_group, dispatch_get_main_queue(), ^{
			m_delegate->context_did_update_pbwt_samples_to_traceback_positions(*this);
		});
	}
	
	
	void segmentation_lp_context::start_update_sample_task(
		std::size_t const lb,
		pbwt_sample &&sample,
		text_position_vector &&right_bounds
	)
	{
		// Start the task.
		// Use pointers to avoid problems if m_update_pbwt_tasks needs to reallocate.
		auto &task_ptr(m_update_pbwt_tasks.emplace_back(new update_pbwt_task(lb, std::move(sample), std::move(right_bounds))));
		auto *task(task_ptr.get());
		dispatch_group_async(*m_update_samples_group, *m_producer_queue, ^{
			task->update_pbwt();
		});
	}
	
	
	void segmentation_lp_context::find_segments_greedy()
	{
		typedef update_pbwt_task::pbwt_sample_vector pbwt_sample_vector;
		
		// Iterate the divergence samples from the second one. Try to join the rightmost segment to the current
		// segment run by calculating its size with respect to the current left bound. If the size is less than the size 
		// calculated with DP, continue. Otherwise, end the current run (to the previously iterated sample) and start a new run.
		std::size_t current_lb(m_update_pbwt_tasks.front()->left_bound()); // Left bound.
		std::size_t prev_size(m_segmentation_traceback_res.front().segment_size);
		auto *prev_sample(&m_update_pbwt_tasks.front()->samples().front());
		
		auto pbwt_samples(m_update_pbwt_tasks	| ranges::view::transform([](auto &task) -> pbwt_sample_vector & { return task->samples(); })
						  						| ranges::view::join);
		for (auto const &tup : ranges::view::zip(pbwt_samples, m_segmentation_traceback_res) | ranges::view::drop(1))
		{
			auto &sample(std::get <0>(tup));
			auto const &dp_arg(std::get <1>(tup));
			assert(sample.rb == dp_arg.rb);
			
			auto const sample_size(sample.context.unique_substring_count_lhs(current_lb));
			if (sample_size <= m_max_segment_size)
				prev_size = sample_size;
			else
			{
				m_reduced_traceback.emplace_back(current_lb, prev_sample->rb, prev_size);
				prev_size = dp_arg.segment_size;
				
				current_lb = prev_sample->rb;
				m_reduced_pbwt_samples.emplace_back(std::move(*prev_sample));
			}
			
			prev_sample = &sample;
		}
		
		m_reduced_traceback.emplace_back(current_lb, prev_sample->rb, prev_size);
		m_reduced_pbwt_samples.emplace_back(std::move(*prev_sample));
		m_update_pbwt_tasks.clear();
	}
	
	
	void segmentation_lp_context::make_cumulative_sum(substring_copy_number_vector &vec) const
	{
		std::size_t cumulative_sum(0);
		for (auto &cn : vec)
		{
			cumulative_sum += cn.copy_number;
			cn.copy_number = cumulative_sum;
		}
	}
	
	
	void segmentation_lp_context::output_gaps(
		std::ostream &ostream,
		std::size_t const text_pos,
		std::size_t const text_length
	) const
	{
		std::fill_n(std::ostreambuf_iterator <char>(ostream), text_length, '-');
	}
	
	
	void segmentation_lp_context::output_substring(
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
	
	
	void segmentation_lp_context::join_segments_and_output(std::ostream &ostream, sequence_vector const &sequences, segment_joining const seg_joining)
	{
		// Count the instances of each substring.
		assert(m_reduced_traceback.size() == m_reduced_pbwt_samples.size());
		substring_copy_number_matrix substring_copy_numbers(m_reduced_traceback.size());

		lb::parallel_for_each(
			ranges::view::zip(m_reduced_traceback, m_reduced_pbwt_samples, substring_copy_numbers),
			[this, seg_joining](auto const &tup){
				auto const &dp_arg(std::get <0>(tup));
				auto const &sample(std::get <1>(tup));
				auto &substring_cn(std::get <2>(tup));
				
				// Count the instances w.r.t. dp_arg’s left bound and sort in decreasing order.
				// Then, in case of non-greedy matching, fill the segment up to the maximum
				// segment size by copying substrings in proportion to their occurrence.
				sample.context.unique_substring_count_idxs_lhs(dp_arg.lb, substring_cn);
				make_cumulative_sum(substring_cn);
			}
		);
		
		// *this may be invalid after calling context_did_output_founders().
		switch (seg_joining)
		{
			case segment_joining::GREEDY:
				join_greedy_and_output(ostream, substring_copy_numbers, sequences);
				m_delegate->context_did_output_founders(*this);
				break;
				
			case segment_joining::BIPARTITE_MATCHING:
				join_with_bipartite_matching_and_output(ostream, substring_copy_numbers, sequences);
				break;
				
			case segment_joining::RANDOM:
				join_random_order_and_output(ostream, substring_copy_numbers, sequences);
				m_delegate->context_did_output_founders(*this);
				break;
				
			case segment_joining::PBWT_ORDER:
				join_pbwt_order_and_output(ostream, substring_copy_numbers, sequences);
				m_delegate->context_did_output_founders(*this);
				break;
				
			default:
				libbio_fail("Unexpected segment joining method.");
				break;
		}
	}
	
	
	std::pair <std::uint32_t, std::uint8_t> segmentation_lp_context::init_permutations(permutation_vector &permutations) const
	{
		// Since we output row-wise, all of the permutation vectors need to be pre-calculated.
		// Try to save some space by using SDSL’s integer vectors.
		// seq_count is not a string index but the following will allow storing it.
		auto const seq_count(m_pbwt_ctx.size());
		auto const bits_needed(lb::bits::highest_bit_set(seq_count));
		auto const permutation_max((1 << bits_needed) - 1);
		assert(seq_count <= permutation_max);
		
		for (auto &permutation : permutations)
		{
			sdsl::int_vector <0> temp_permutation(m_max_segment_size, 0, bits_needed);
			permutation = std::move(temp_permutation);
		}
		
		return std::make_pair(permutation_max, bits_needed);
	}
	
	
	// Fill index_pairs, to_lhs_substring, to_rhs_string.
	std::pair <std::size_t, std::size_t>
	segmentation_lp_context::greedy_create_index_pairs(
		substring_copy_number_vector const &lhs_cn_vector,
		substring_copy_number_vector const &rhs_cn_vector,
		pbwt_sample const &lhs_sample,
		pbwt_sample const &rhs_sample,
		sdsl::int_vector <0> const &rhs_matching,
		std::uint64_t const permutation_max,
		std::vector <detail::substring_index_pair> &index_pairs,
		std::vector <detail::substring_index_pair> &index_pairs_buffer,
		sdsl::int_vector <0> &to_lhs_substring,
		sdsl::int_vector <0> &to_rhs_string
	) const
	{
		std::size_t lhs_substring_count(0);
		std::size_t rhs_substring_count(0);
		
		// Fill the pairs of matched string indices, count the occurrence of each pair and sort by it.
		// Also update to_lhs_substring s.t. string numbers can be converted to unique substring numbers
		// and to_rhs_string s.t. unique substring numbers can be converted to string numbers.
		assert(0 < rhs_matching.size());
		auto const seq_count(m_pbwt_ctx.size());
		index_pairs.clear();
		index_pairs.resize(seq_count);	// Fill with the default constructor.
		auto const &lhs_input_permutation(lhs_sample.context.input_permutation());
		auto const &rhs_input_permutation(rhs_sample.context.input_permutation());
		auto lhs_cn_it(lhs_cn_vector.begin());
		auto rhs_cn_it(rhs_cn_vector.begin());
		std::size_t copy_number(0);
		std::size_t lhs_substring_number(0);
		std::size_t rhs_substring_number(0);
		std::size_t lhs_substring_number_conv(rhs_matching[0]);
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
				index_pairs[lhs_string_number].lhs_idx = lhs_substring_number_conv;
				index_pairs[rhs_string_number].rhs_idx = rhs_substring_number;
				
				// Update to_lhs_substring and to to_rhs_string.
				// Note: to_rhs_string gets maximum indices, not minimum. At this time it doesn’t matter, though.
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
		std::fill(to_rhs_string.begin() + rhs_substring_count, to_rhs_string.end(), permutation_max);
		
		assert(lhs_cn_vector.end() == ++lhs_cn_it);
		assert(rhs_cn_vector.end() == ++rhs_cn_it);
		
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
	void segmentation_lp_context::greedy_create_matching(
		std::vector <detail::substring_index_pair> const &index_pairs,
		std::uint64_t const matching_max,
		sdsl::int_vector <0> &lhs_matching,
		sdsl::int_vector <0> &rhs_matching
	) const
	{
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
			while (lhs_idx < m_max_segment_size && lhs_matching[lhs_idx] != matching_max)
				++lhs_idx;
			
			while (rhs_idx < m_max_segment_size && rhs_matching[rhs_idx] != matching_max)
				++rhs_idx;
			
			if (lhs_idx == m_max_segment_size || rhs_idx == m_max_segment_size)
			{
				assert(lhs_idx == rhs_idx);
				break;
			}
			
			lhs_matching[lhs_idx] = rhs_idx;
			rhs_matching[rhs_idx] = lhs_idx;
		}
	}

	
	void segmentation_lp_context::join_greedy_and_output(
		std::ostream &stream,
		substring_copy_number_matrix const &substrings_to_output,
		sequence_vector const &sequences
	) const
	{
		// Use the greedy algorithm to generate the permutations.
		permutation_vector permutations(m_reduced_traceback.size());
		auto const res(init_permutations(permutations));
		auto const permutation_max(std::get <0>(res));
		auto const permutation_bits_needed(std::get <1>(res));
		assert(m_max_segment_size);
		auto const matching_bits_needed(lb::bits::highest_bit_set(m_max_segment_size - 1));
		auto const matching_max((1 << matching_bits_needed) - 1);
		auto const seq_count(m_pbwt_ctx.size());

		std::vector <detail::substring_index_pair> index_pairs;			// Pairs of matched string indices.
		std::vector <detail::substring_index_pair> index_pairs_buffer;	// For radix sorting and eliminating duplicates.
		sdsl::int_vector <0> lhs_matching(m_max_segment_size, matching_max, matching_bits_needed);
		sdsl::int_vector <0> rhs_matching(m_max_segment_size, 0, matching_bits_needed);
		sdsl::int_vector <0> to_lhs_substring(seq_count, 0, matching_bits_needed);		// For converting the string numbers to unique substring numbers.
		sdsl::int_vector <0> to_rhs_string(1 + matching_max, 0, permutation_bits_needed);

		// Set the initial state with the substring numbers and the minimum unique substring indices.
		// rhs_matching is used first.
		std::iota(rhs_matching.begin(), rhs_matching.end(), 0);	// Start with identity.
		{
			auto const &cn_vec(substrings_to_output.front());
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
				permutation[i] = permutation_max;
				++i;
			}
		}
		
		for (auto const &tup : ranges::view::zip(substrings_to_output, m_reduced_pbwt_samples, permutations) | ranges::view::sliding(2))
		{
			auto const &lhs(tup[0]);
			auto const &rhs(tup[1]);
			auto const &lhs_cn_vector(std::get <0>(lhs));
			auto const &rhs_cn_vector(std::get <0>(rhs));
			auto const &lhs_sample(std::get <1>(lhs));
			auto const &rhs_sample(std::get <1>(rhs));
			auto const &lhs_output_permutation(std::get <2>(lhs));
			auto &rhs_output_permutation(std::get <2>(rhs));
			
			// Create pairs of string indices of the unique substrings and sort them.
			auto const res(greedy_create_index_pairs(
				lhs_cn_vector,
				rhs_cn_vector,
				lhs_sample,
				rhs_sample,
				rhs_matching,
				permutation_max,
				index_pairs,
				index_pairs_buffer,
				to_lhs_substring,
				to_rhs_string
			));
			std::size_t lhs_substring_count(res.first);
			std::size_t rhs_substring_count(res.second);
			greedy_create_matching(index_pairs, matching_max, lhs_matching, rhs_matching);
			
			{
				// Create the next permutation.
				std::size_t idx(0);
				std::size_t next_gap_substring_idx(lhs_substring_count);
				for (auto const lhs_string_idx : lhs_output_permutation)
				{
					// Check if there is a gap sequence on the left.
					// Join to the gap sequences in an arbitrary order.
					std::size_t lhs_substring_idx(0);
					if (lhs_string_idx == permutation_max)
						lhs_substring_idx = next_gap_substring_idx++;
					else
						lhs_substring_idx = to_lhs_substring[lhs_string_idx];
					
					auto const matched_idx(lhs_matching[lhs_substring_idx]);
					auto const rhs_string_idx(to_rhs_string[matched_idx]);
					rhs_output_permutation[idx++] = rhs_string_idx;
				}
			}
		}
		
		output_in_permutation_order(stream, permutations, sequences, permutation_max);
	}
	
	
	void segmentation_lp_context::join_with_bipartite_matching_and_output(
		std::ostream &stream,
		substring_copy_number_matrix const &substrings_to_output,
		sequence_vector const &sequences
	)
	{
		// FIXME: write me.
		m_delegate->context_did_output_founders(*this);
	}
	
	
	void segmentation_lp_context::join_random_order_and_output(
		std::ostream &stream,
		substring_copy_number_matrix const &substrings_to_output,
		sequence_vector const &sequences
	) const
	{
		permutation_vector permutations(m_reduced_traceback.size());
		auto const res(init_permutations(permutations));

		// Instantiate an std::mt19937 with the given seed and generate a permutation of the indices for each segment.
		std::mt19937 urbg(m_random_seed);
		
		// Fill the permutations and shuffle.
		for (auto const &tup : ranges::view::zip(substrings_to_output, permutations))
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
		
		output_in_permutation_order(stream, permutations, sequences, res.first);
	}
	
	
	void segmentation_lp_context::join_pbwt_order_and_output(
		std::ostream &stream,
		substring_copy_number_matrix const &substrings_to_output,
		sequence_vector const &sequences
	) const
	{
		// Create iterators to each vector of runs.
		std::vector <substring_copy_number_vector::const_iterator> cn_iterators(substrings_to_output.size());
		for (auto const &tup : ranges::view::zip(substrings_to_output, cn_iterators))
			std::get <1>(tup) = std::get <0>(tup).cbegin();
		
		// Output.
		for (std::size_t row(0); row < m_max_segment_size; ++row)
		{
			for (auto const &tup : ranges::view::zip(cn_iterators, m_reduced_traceback))
			{
				auto &it(std::get <0>(tup));
				auto const &dp_arg(std::get <1>(tup));
				auto const ending_index(it->copy_number); // Running sum.
				
				if (row == ending_index)
					++it;
				
				assert(row < it->copy_number);
				
				auto const substring_idx(it->substring_idx);
				output_substring(stream, substring_idx, dp_arg.lb, dp_arg.text_length(), sequences);
			}
			stream << '\n';
		}
		
#ifndef NDEBUG
		for (auto const &tup : ranges::view::zip(substrings_to_output, cn_iterators))
			assert(++std::get <1>(tup) == std::get <0>(tup).cend());
#endif
	}
	
	
	void segmentation_lp_context::output_in_permutation_order(
		std::ostream &stream,
		permutation_vector const &permutations,
		sequence_vector const &sequences,
		std::uint32_t const permutation_max
	) const
	{
		// Output.
		for (std::size_t row(0); row < m_max_segment_size; ++row)
		{
			for (auto const &tup : ranges::view::zip(permutations, m_reduced_traceback))
			{
				auto &permutation(std::get <0>(tup));
				auto const &dp_arg(std::get <1>(tup));
				auto const substring_idx(permutation[row]);
				if (permutation_max == substring_idx)
					output_gaps(stream, dp_arg.lb, dp_arg.text_length());
				else
					output_substring(stream, substring_idx, dp_arg.lb, dp_arg.text_length(), sequences);
			}
			stream << '\n';
		}
	}
	
	
	void segmentation_lp_context::output_segmentation_status(std::size_t const j) const
	{
		if (0 == j % 100000)
			std::cerr << ' ' << j << std::flush;
	}
	
	
	void segmentation_lp_context::output_segmentation_status_mq(std::size_t const j) const
	{
		if (0 == j % 100000)
		{
			dispatch_async(dispatch_get_main_queue(), ^{
				std::cerr << ' ' << j << std::flush;
			});
		}
	}
	
	
	void calculate_segmentation_lp_dp_arg(
		typename buffering_pbwt_context::divergence_count_list const &divergence_value_counts,
		segmentation_traceback_vector const &segmentation_traceback_dp,
		segmentation_traceback_vector_rmq const &segmentation_traceback_dp_rmq,
		std::size_t const seq_count,
		std::size_t const segment_length,
		std::size_t const lb,	// Inclusive.
		std::size_t const text_pos,
		segmentation_dp_arg &min_arg
	)
	{
		// Take the current divergence values and their counts in divergence value order.
		auto it(divergence_value_counts.cbegin_pairs());
		auto const end(divergence_value_counts.cend_pairs());
		assert(it != end);
		
		// Find the minimum, similar to Ukkonen's equation 1.
		std::size_t segment_size_diff(it->second);
		
		std::size_t dp_lb(lb);
		std::size_t dp_rb(it->first); // Initial value not used.
		
		// Consider the whole range in case the segment size is smaller than seq_count.
		if (lb == dp_rb)
		{
			auto const segment_size(seq_count - segment_size_diff);
			segmentation_dp_arg const current_arg(lb, 1 + text_pos, segment_size, segment_size);
			if (current_arg < min_arg)
				min_arg = current_arg;
			
			++it;
			dp_rb = it->first;
			segment_size_diff += it->second;
		}
		
		++it;
		while (true)
		{
			if (it == end)
				break;
			
			// Range of possible cutting points.
			dp_lb = dp_rb;
			dp_rb = it->first;
			std::size_t dp_rb_c(dp_rb);
			assert(0 != it->first);
			assert(dp_lb < dp_rb);
			assert(dp_rb <= 1 + text_pos);
			
			// Check if the range needs to and can be contracted.
			// First, verify that the current segment fits after the DP search range.
			if (text_pos + 2 - segment_length < dp_rb_c)
				dp_rb_c = text_pos + 2 - segment_length;
			
			// Second, verify that one segment fits before the DP search range.
			// (Considering segment end positions; it must hold that lb + segment_length ≤ dp_lb.)
			if (dp_lb < lb + segment_length)
			{
				if (lb + segment_length < dp_rb_c)
					dp_lb = lb + segment_length;
				else
					goto continue_loop;
			}
			
			// Check that the range is still valid.
			if (dp_lb < dp_rb_c)
			{
				// Convert to segmentation_traceback indexing.
				assert(segment_length <= dp_lb);
				assert(segment_length <= dp_rb_c);
				auto const dp_lb_tb(dp_lb - segment_length);
				auto const dp_rb_tb(dp_rb_c - segment_length);
				auto const idx(segmentation_traceback_dp_rmq(dp_lb_tb, dp_rb_tb));
				
				auto const &boundary_segment(segmentation_traceback_dp[idx]);
				auto const lhs(boundary_segment.segment_max_size);
				auto const rhs(seq_count - segment_size_diff);
				
				segmentation_dp_arg current_arg(idx + segment_length, 1 + text_pos, lb::max_ct(lhs, rhs), rhs);
				if (current_arg < min_arg)
					min_arg = current_arg;
			}
			
			// Next iteration.
		continue_loop:
			segment_size_diff += it->second;
			++it;
		}
	}
}
