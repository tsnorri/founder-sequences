/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <founder_sequences/segmentation_lp_context.hh>
#include <libbio/algorithm.hh>

namespace lb = libbio;


namespace {
	static uint32_t const s_output_status_freq = 1000;
}


namespace founder_sequences {
	
	void calculate_segmentation_lp_dp_arg(
		typename buffering_pbwt_context::divergence_count_list const &divergence_value_counts,
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
		dispatch_async(*m_producer_queue, ^{
			m_pbwt_ctx.set_sample_rate(m_delegate->pbwt_sample_rate());
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
			
			m_pbwt_ctx.process <false>(
				segment_length - 1,
				[this](std::size_t const idx){ output_segmentation_status_mq(1 + idx, m_pbwt_ctx.samples().size()); },
				[this, lb, rb](){ generate_traceback_part_2(lb, rb); },
				[this](void (^block)()){ m_dispatch_helper->dispatch(*m_producer_queue, block); }
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
			
			m_pbwt_ctx.process <true>(
				limit,
				[this, lb, seq_count, segment_length](std::size_t const idx, std::size_t const buffer_idx, dispatch_semaphore_t process_sema){
					// Counts is just a column object, copy.
					auto const &counts(m_pbwt_ctx.divergence_value_counts(buffer_idx));
					auto const sample_count(m_pbwt_ctx.samples().size());
					m_dispatch_helper->dispatch(*m_consumer_queue, ^{
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

						dispatch_semaphore_signal(process_sema);
						output_segmentation_status_mq(1 + idx, sample_count);
					});
				},
				[this, lb, rb](){ generate_traceback_part_3(lb, rb); },
				[this](void (^block)()){ m_dispatch_helper->dispatch(*m_producer_queue, block); }
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
			
			m_pbwt_ctx.process <true>(
				limit,
				[this, lb, seq_count, segment_length](std::size_t const idx, std::size_t const buffer_idx, dispatch_semaphore_t process_sema){
					auto const &counts(m_pbwt_ctx.divergence_value_counts(buffer_idx));
					auto const sample_count(m_pbwt_ctx.samples().size());
					m_dispatch_helper->dispatch(*m_consumer_queue, ^{
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
						
						dispatch_semaphore_signal(process_sema);
						output_segmentation_status_mq(1 + idx, sample_count);
					});
				},
				[this, lb, rb](){ generate_traceback_part_4(lb, rb); },
				[this](void (^block)()){ m_dispatch_helper->dispatch(*m_producer_queue, block); }
			);
		});
	}
	
	
	void segmentation_lp_context::generate_traceback_part_4(std::size_t const lb, std::size_t const rb)
	{
		// Calculate the size of the final segment.
		dispatch_async(*m_producer_queue, ^{
			auto const seq_count(m_pbwt_ctx.size());
			auto const segment_length(m_delegate->segment_length());
			
			m_pbwt_ctx.process <false>(
				rb,
				[this](std::size_t const idx){ output_segmentation_status_mq(1 + idx, m_pbwt_ctx.samples().size()); },
				[this, lb, rb, seq_count, segment_length](){
					m_dispatch_helper->dispatch(dispatch_get_main_queue(), ^{
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
				},
				[this](void (^block)()){ m_dispatch_helper->dispatch(*m_producer_queue, block); }
			);
		});
	}
	
	
	void segmentation_lp_context::follow_traceback()
	{
		// Follow the traceback.
		m_dispatch_helper->dispatch(dispatch_get_main_queue(), ^{
			lb::log_time(std::cerr);
			std::cerr << "Following the traceback…" << std::flush;
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
		
		auto const segment_count(m_segmentation_traceback_res.size());
		dispatch_async(dispatch_get_main_queue(), ^{
			std::cerr << " there were " << segment_count << " segments the maximum size of which was " << m_max_segment_size << '.' << std::endl;
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
				assert(last_sample.rb <= right_bounds.front());
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
		pbwt_sample_type &&sample,
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
	
	
	void segmentation_lp_context::find_segments_greedy(segmentation_container &container)
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
				container.reduced_traceback.emplace_back(current_lb, prev_sample->rb, prev_size);
				prev_size = dp_arg.segment_size;
				
				current_lb = prev_sample->rb;
				container.reduced_pbwt_samples.emplace_back(std::move(*prev_sample));
			}
			
			prev_sample = &sample;
		}
		
		container.max_segment_size = m_max_segment_size;
		container.reduced_traceback.emplace_back(current_lb, prev_sample->rb, prev_size);
		container.reduced_pbwt_samples.emplace_back(std::move(*prev_sample));
		m_update_pbwt_tasks.clear();
	}
	
	
	void segmentation_lp_context::output_segmentation_status(std::size_t const j, std::size_t const sample_count) const
	{
		if (0 == j % s_output_status_freq)
			output_segmentation_status_2(j, sample_count);
	}
	
	
	void segmentation_lp_context::output_segmentation_status_mq(std::size_t const j, std::size_t const sample_count) const
	{
		if (0 == j % s_output_status_freq)
		{
			m_dispatch_helper->dispatch(dispatch_get_main_queue(), ^{
				output_segmentation_status_2(j, sample_count);
			});
		}
	}
	
	
	void segmentation_lp_context::output_segmentation_status_2(std::size_t const j, std::size_t const sample_count) const
	{
		std::cerr << ' ' << j << ", " << sample_count;
		if (1 == sample_count)
			std::cerr << " sample;";
		else
			std::cerr << " samples;";
		std::cerr << std::flush;
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
