/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <founder_sequences/bipartite_matcher.hh>
#include <founder_sequences/create_segment_texts_task.hh>
#include <founder_sequences/merge_segments_task.hh>
#include <founder_sequences/segment_text.hh>


namespace lb	= libbio;


namespace founder_sequences
{
	void bipartite_matcher::match()
	{
		lb::dispatch_ptr <dispatch_group_t> group(dispatch_group_create());
		
		auto const seq_count(m_delegate->sequence_count());
		auto const max_segment_size(m_delegate->max_segment_size());
		auto const &pbwt_samples(m_delegate->pbwt_samples());
		
		// Start tasks for creating segment_texts.
		m_segment_texts.resize(pbwt_samples.size());
		for (auto const &tup : ranges::view::zip(pbwt_samples, *m_substrings_to_output, m_segment_texts))
		{
			auto const &sample(std::get <0>(tup));
			auto const &copy_numbers(std::get <1>(tup));
			auto &segment_texts(std::get <2>(tup));
			
			auto &task_ptr(m_tasks.emplace_back(new create_segment_texts_task(sample, copy_numbers, max_segment_size, seq_count, segment_texts)));
			// The dispatch group is retained.
			// Use a raw pointer in case m_tasks needs to reallocate instead of using a reference to the inserted unique_ptr.
			// https://clang.llvm.org/docs/BlockLanguageSpec.html#c-extensions
			auto *task(task_ptr.get());
			assert(task);
			dispatch_group_async(*group, *m_producer_queue, ^{
				assert(task);
				task->execute();
			});
		}
		
		dispatch_group_notify(*group, *m_consumer_queue, ^{
			start_matching_tasks();
		});
	}
	
	
	void bipartite_matcher::output_segments(std::ostream &stream, sequence_vector const &sequences)
	{
		::founder_sequences::output_segments(
			stream,
			m_delegate->reduced_traceback(),
			m_segment_texts,
			sequences
		);
	}
	
	
	void bipartite_matcher::start_matching_tasks()
	{
		// All the previous tasks have now finished executing.
		m_tasks.clear();
		
		// Create the initial permutation.
		create_initial_permutation();
		
		auto const task_count(m_segment_texts.size() - 1);
		
		m_matchings.clear();
		m_task_completion.clear();
		m_started_from.clear();
		
		m_matchings.resize(task_count);
		
		{
			using std::swap;
			decltype(m_task_completion) temp(task_count);
			swap(m_task_completion, temp);
		}
		{
			using std::swap;
			decltype(m_started_from) temp(task_count);
			swap(m_started_from, temp);
		}
		
		// Create the merging tasks.
		auto const set_scoring_method(m_delegate->bipartite_set_scoring_method());
		lb::dispatch_ptr <dispatch_group_t> group(dispatch_group_create());
		std::size_t task_idx(0);
		for (auto const &pair : m_segment_texts | ranges::view::sliding(2))
		{
			auto &matching(m_matchings[task_idx]);
			auto &task_ptr(m_tasks.emplace_back(new merge_segments_task(task_idx++, *this, pair[0], pair[1], matching, set_scoring_method)));
			auto *task(task_ptr.get());
			dispatch_group_async(*group, *m_producer_queue, ^{
				task->execute();
			});
		}
		assert(task_count == m_tasks.size());
		
		dispatch_group_notify(*group, *m_consumer_queue, ^{
			m_delegate->matcher_did_finish(*this);
		});
	}
	
	
	void bipartite_matcher::task_did_finish(merge_segments_task &task)
	{
		auto idx(task.task_index());
		assert(idx < m_task_completion.size());
		m_task_completion[idx].fetch_or(1);
		if (idx == m_expecting)
		{
			while (0 == m_started_from[idx].fetch_or(1))
			{
				// Critical section.
				auto const next_to_be_handled(create_permutations(idx));
				// End of critical section.
				
				// Check if we’re done.
				if (next_to_be_handled == m_task_completion.size())
					return;
				
				// Set the next starting point.
				// After this, other threads may enter this loop.
				auto const prev_expected(m_expecting.exchange(next_to_be_handled));
				assert(prev_expected < next_to_be_handled);
			
				// See if the task at the next starting point was completed meanwhile.
				if (!m_task_completion[next_to_be_handled])
					return;
			
				// The task at the starting point was actually completed.
				// See if permutation creation was already started by running the loop again.
				idx = next_to_be_handled;
			}
		}
	}
	
	
	void bipartite_matcher::create_initial_permutation()
	{
		// Store the segment text order.
		m_segment_text_permutation.clear();
		m_segment_text_permutation.resize(m_delegate->max_segment_size());
		std::iota(m_segment_text_permutation.begin(), m_segment_text_permutation.end(), 0);
		
		// Create the string permutation.
		auto const &segment_texts(m_segment_texts.front());
		auto &permutations(m_delegate->permutations());
		auto &first_permutation(permutations.front());
		assert(segment_texts.size() == first_permutation.size());
		
		// For each segment_text, store a representative into first_permutation.
		std::size_t i(0);
		for (auto const &seg_text : segment_texts)
		{
			auto const &non_copied_seg_text(segment_texts[seg_text.row_number(i)]);
			assert(!non_copied_seg_text.is_copied());
			first_permutation[i] = non_copied_seg_text.first_sequence_index();
			
			++i;
		}
	}
	
	
	std::size_t bipartite_matcher::create_permutations(std::size_t idx)
	{
		auto &permutations(m_delegate->permutations());
		while (m_task_completion[idx])
		{
			auto const &matching(m_matchings[idx]);
			auto const &segment_texts(m_segment_texts[1 + idx]);
			auto &permutation(permutations[1 + idx]);
			assert(permutation.size() == m_segment_text_permutation.size());
			
			// Fill the next permutation and update m_segment_text_permutation.
			std::size_t i(0);
			for (auto &seg_idx : m_segment_text_permutation)
			{
				assert(i < permutation.size());
				
				auto const matched_seg_idx(matching[seg_idx]);
				auto const &seg_text(segment_texts[matched_seg_idx]);
				auto const &non_copied_seg_text(segment_texts[seg_text.row_number(i)]);
				assert(!non_copied_seg_text.is_copied());
				permutation[i] = non_copied_seg_text.first_sequence_index();
				
				seg_idx = matched_seg_idx;
				++i;
			}
			
			++idx;
		}
		
		return idx;
	}
}
