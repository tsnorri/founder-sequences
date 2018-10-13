/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <founder_sequences/segmentation_lp_context.hh>
#include <libbio/algorithm.hh>

namespace lb = libbio;


namespace founder_sequences {
	
	void calculate_segmentation_lp_dp_arg(
		typename pbwt_context::divergence_count_list const &divergence_value_counts,
		segmentation_traceback_vector const &segmentation_traceback_dp,
		segmentation_traceback_vector_rmq const &segmentation_traceback_dp_rmq,
		std::size_t const seq_count,
		std::size_t const segment_length,
		std::size_t const lb,	// Inclusive.
		std::size_t const text_pos,
		segmentation_dp_arg &min_arg
	);
	
	
	void segmentation_lp_context::generate_traceback(std::size_t const lb, std::size_t const rb)
	{
		// Calculate the first L - 1 columns, which gives the required result for calculating M(L).
		// idx is 0-based, m_segment_length is 1-based.
		m_step_max = m_delegate->sequences().front().size();
		
		dispatch_async(*m_producer_queue, ^{
			m_pbwt_ctx.set_sample_rate(m_delegate->pbwt_sample_rate());
			m_pbwt_ctx.prepare();
			
			auto const seq_length(m_pbwt_ctx.sequence_length());
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
			
			m_pbwt_ctx.process <lb::pbwt::context_field::DIVERGENCE_VALUE_COUNTS>(
				segment_length - 1,
				[this](){
					m_current_step.store(m_pbwt_ctx.sequence_idx() + 1, std::memory_order_relaxed);
					m_current_pbwt_sample_count.store(m_pbwt_ctx.samples().size(), std::memory_order_relaxed);
				}
			);
			
			generate_traceback_part_2(lb, rb);
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
			
			m_pbwt_ctx.process <lb::pbwt::context_field::DIVERGENCE_VALUE_COUNTS>(
				limit,
				[this, lb, seq_count, segment_length](){
					auto const idx(m_pbwt_ctx.sequence_idx());
					auto const &counts(m_pbwt_ctx.output_divergence_value_counts());
					auto const sample_count(m_pbwt_ctx.samples().size());
					
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
					
					m_current_step.store(1 + idx, std::memory_order_relaxed);
					m_current_pbwt_sample_count.store(sample_count, std::memory_order_relaxed);
				}
			);
			
			generate_traceback_part_3(lb, rb);
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
			
			m_pbwt_ctx.process <lb::pbwt::context_field::DIVERGENCE_VALUE_COUNTS>(
				limit,
				[this, lb, seq_count, segment_length](){
					auto const idx(m_pbwt_ctx.sequence_idx());
					auto const &counts(m_pbwt_ctx.output_divergence_value_counts());
					auto const sample_count(m_pbwt_ctx.samples().size());
					
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
					
					m_current_step.store(1 + idx, std::memory_order_relaxed);
					m_current_pbwt_sample_count.store(sample_count, std::memory_order_relaxed);
				}
			);
			
			generate_traceback_part_4(lb, rb);
		});
	}
	
	
	void segmentation_lp_context::generate_traceback_part_4(std::size_t const lb, std::size_t const rb)
	{
		// Calculate the size of the final segment.
		dispatch_async(*m_producer_queue, ^{
			auto const seq_count(m_pbwt_ctx.size());
			auto const segment_length(m_delegate->segment_length());
			
			m_pbwt_ctx.process <lb::pbwt::context_field::DIVERGENCE_VALUE_COUNTS>(
				rb,
				[this](){
					auto const idx(m_pbwt_ctx.sequence_idx());
					m_current_step.store(1 + idx, std::memory_order_relaxed);
					m_current_pbwt_sample_count.store(m_pbwt_ctx.samples().size(), std::memory_order_relaxed);
				}
			);
			
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
			m_current_step.store(m_step_max, std::memory_order_relaxed);
			
			follow_traceback();
		});
	}
	
	
	void segmentation_lp_context::follow_traceback()
	{
		// Follow the traceback.
		m_delegate->context_will_follow_traceback(*this);
		
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
		
		auto const segment_count(m_segmentation_traceback_res.size());
		dispatch_async(dispatch_get_main_queue(), ^{
			m_delegate->context_did_finish_traceback(*this, segment_count, m_max_segment_size);
		});
	}
	
	
	// Given a collection of PBWT samples, update them s.t. their right bound matches that of
	// the traceback arguments.
	void segmentation_lp_context::update_samples_to_traceback_positions()
	{
		dispatch_async(*m_producer_queue, ^{
			m_update_samples_group.reset(dispatch_group_create());
			
			auto &pbwt_samples(m_pbwt_ctx.samples());
			auto const sample_count(pbwt_samples.size());
			auto traceback_it(m_segmentation_traceback_res.cbegin());
			auto const traceback_end(m_segmentation_traceback_res.cend());
			
			m_current_step = 0;
			m_step_max = sample_count;
			
			m_delegate->context_will_start_update_samples_tasks(*this);
			
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
				while (traceback_it < traceback_end && traceback_it->rb < next_sample.sequence_idx())
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
				else
				{
					// Did not start a task.
					m_current_step.fetch_add(1, std::memory_order_relaxed);
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
					assert(last_sample.sequence_idx() <= right_bounds.front());
					start_update_sample_task(lb, std::move(last_sample), std::move(right_bounds));
				}
				else
				{
					m_current_step.fetch_add(1, std::memory_order_relaxed);
				}
			}
			
			m_delegate->context_did_start_update_samples_tasks(*this);
			
			// Remove the remaining samples.
			pbwt_samples.clear();
			
			dispatch_group_notify(*m_update_samples_group, dispatch_get_main_queue(), ^{
				m_delegate->context_did_update_pbwt_samples_to_traceback_positions(*this);
			});
		});
	}
	
	
	void segmentation_lp_context::start_update_sample_task(
		std::size_t const lb,
		pbwt_sample_type &&sample,
		text_position_vector &&right_bounds
	)
	{
		// Start the task.
		// Use pointers to avoid problems if m_update_pbwt_tasks needs to reallocate.
		auto &task_ptr(m_update_pbwt_tasks.emplace_back(new update_pbwt_task(*this, lb, std::move(sample), std::move(right_bounds))));
		auto *task(task_ptr.get());
		dispatch_group_async(*m_update_samples_group, *m_producer_queue, ^{
			task->execute();
		});
	}
	
	
	void segmentation_lp_context::find_segments_greedy()
	{
		dispatch_async(*m_producer_queue, ^{
			
			typedef update_pbwt_task::pbwt_sample_vector pbwt_sample_vector;
			
			segmentation_container container;
			
			// Iterate the divergence samples from the second one. Try to join the rightmost segment to the current
			// segment run by calculating its size with respect to the current left bound. If the size is less than the size 
			// calculated with DP, continue. Otherwise, end the current run (to the previously iterated sample) and start a new run.
			std::size_t current_lb(m_update_pbwt_tasks.front()->left_bound()); // Left bound.
			std::size_t prev_size(m_segmentation_traceback_res.front().segment_size);
			auto *prev_sample(&m_update_pbwt_tasks.front()->samples().front());
			
			// Progress tracking.
			m_current_step = 0;
			m_step_max = m_segmentation_traceback_res.size();
			m_delegate->context_will_merge_segments(*this);
			
			auto pbwt_samples(m_update_pbwt_tasks	| ranges::view::transform([](auto &task) -> pbwt_sample_vector & { return task->samples(); })
							  						| ranges::view::join);
			for (auto const &tup : ranges::view::zip(pbwt_samples, m_segmentation_traceback_res) | ranges::view::drop(1))
			{
				auto &sample(std::get <0>(tup));
				auto const &dp_arg(std::get <1>(tup));
				assert(sample.sequence_idx() == dp_arg.rb);
				
				auto const sample_size(sample.unique_substring_count_lhs(current_lb));
				if (sample_size <= m_max_segment_size)
					prev_size = sample_size;
				else
				{
					container.reduced_traceback.emplace_back(current_lb, prev_sample->sequence_idx(), prev_size);
					prev_size = dp_arg.segment_size;
					
					current_lb = prev_sample->sequence_idx();
					container.reduced_pbwt_samples.emplace_back(std::move(*prev_sample));
				}
				
				prev_sample = &sample;
				m_current_step.fetch_add(1, std::memory_order_relaxed);
			}
			
			container.max_segment_size = m_max_segment_size;
			container.reduced_traceback.emplace_back(current_lb, prev_sample->sequence_idx(), prev_size);
			container.reduced_pbwt_samples.emplace_back(std::move(*prev_sample));
			m_update_pbwt_tasks.clear();
			
			m_current_step.fetch_add(1, std::memory_order_relaxed);
			
			lb::dispatch_async_fn(dispatch_get_main_queue(), [this, container{std::move(container)}]() mutable {
				m_delegate->context_did_merge_segments(*this, std::move(container));
			});
		});
	}
	
	
	void calculate_segmentation_lp_dp_arg(
		typename pbwt_context::divergence_count_list const &divergence_value_counts,
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
