/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef FOUNDER_SEQUENCES_SEGMENTATION_LP_CONTEXT_HH
#define FOUNDER_SEQUENCES_SEGMENTATION_LP_CONTEXT_HH

#include <founder_sequences/bipartite_matcher.hh>
#include <founder_sequences/founder_sequences.hh>
#include <founder_sequences/greedy_matcher.hh>
#include <founder_sequences/segmentation_container.hh>
#include <founder_sequences/segmentation_dp_arg.hh>
#include <founder_sequences/substring_copy_number.hh>
#include <founder_sequences/update_pbwt_task.hh>
#include <libbio/dispatch_fn.hh>


namespace founder_sequences { namespace detail {
	struct dispatch_helper
	{
		virtual ~dispatch_helper() {}
		virtual void dispatch(dispatch_queue_t queue, void (^block)()) = 0;
	};
	
	struct dispatch_helper_st final : public dispatch_helper
	{
		void dispatch(dispatch_queue_t queue, void (^block)()) override { block(); }
	};
	
	struct dispatch_helper_mt final : public dispatch_helper
	{
		void dispatch(dispatch_queue_t queue, void (^block)()) override { dispatch_async(queue, block); }
	};
}}


namespace founder_sequences {
	
	class segmentation_lp_context;
	
	
	struct segmentation_lp_context_delegate : public virtual segmentation_context_delegate
	{
		virtual std::size_t segment_length() const = 0;
		virtual std::uint64_t pbwt_sample_rate() const = 0;
		virtual alphabet_type const &alphabet() const = 0;
		virtual sequence_vector const &sequences() const = 0;
		virtual void context_did_finish_traceback(segmentation_lp_context &ctx) = 0;
		virtual void context_did_update_pbwt_samples_to_traceback_positions(segmentation_lp_context &ctx) = 0;
	};
	
	
	class segmentation_lp_context final : public segmentation_context
	{
	protected:
		typedef std::vector <std::size_t>					text_position_vector;

	protected:
		buffering_pbwt_context								m_pbwt_ctx;
		segmentation_traceback_vector						m_segmentation_traceback_res;
		segmentation_traceback_vector						m_segmentation_traceback_dp;
		segmentation_traceback_vector_rmq					m_segmentation_traceback_dp_rmq;
		std::uint32_t										m_max_segment_size{};
		
		libbio::dispatch_ptr <dispatch_queue_t>				m_producer_queue;	// May be parallel.
		libbio::dispatch_ptr <dispatch_queue_t>				m_consumer_queue;	// Needs to be serial.
		
		// For updating the PBWT samples.
		std::vector <std::unique_ptr <update_pbwt_task>>	m_update_pbwt_tasks;
		libbio::dispatch_ptr <dispatch_group_t>				m_update_samples_group;
		
		// For matching.
		substring_copy_number_matrix						m_substring_copy_numbers;
		std::unique_ptr <matcher>							m_matcher;
		std::uint32_t										m_permutation_max{};
		std::uint8_t										m_permutation_bits_needed{};
		
		std::unique_ptr <detail::dispatch_helper>			m_dispatch_helper;
		segmentation_lp_context_delegate					*m_delegate{};
		
	public:
		segmentation_lp_context(
			segmentation_lp_context_delegate &delegate,
			libbio::dispatch_ptr <dispatch_queue_t> &producer_queue,
			libbio::dispatch_ptr <dispatch_queue_t> &consumer_queue
		):
			m_pbwt_ctx(8, delegate.sequences(), delegate.alphabet()),
			m_producer_queue(producer_queue),
			m_consumer_queue(consumer_queue),
			m_dispatch_helper(
				delegate.should_run_single_threaded() ?
				static_cast <detail::dispatch_helper *>(new detail::dispatch_helper_st()) :
				static_cast <detail::dispatch_helper *>(new detail::dispatch_helper_mt())
			),
			m_delegate(&delegate)
		{
		}
		
		segmentation_traceback_vector &segmentation_traceback() { return m_segmentation_traceback_res; }
		std::uint32_t max_segment_size() const override { return m_max_segment_size; }
		
		void cleanup() { delete this; }
		
		void generate_traceback(std::size_t lb, std::size_t rb);
		void update_samples_to_traceback_positions();
		void find_segments_greedy(segmentation_container &container);
		
	protected:
		void generate_traceback_part_2(std::size_t const lb, std::size_t const rb);
		void generate_traceback_part_3(std::size_t const lb, std::size_t const rb);
		void generate_traceback_part_4(std::size_t const lb, std::size_t const rb);
		
		void follow_traceback();
		
		void start_update_sample_task(
			std::size_t const lb,
			pbwt_sample_type &&sample,
			text_position_vector &&right_bounds
		);
		
		inline void output_segmentation_status(std::size_t const j, std::size_t const sample_count) const;
		inline void output_segmentation_status_mq(std::size_t const j, std::size_t const sample_count) const;
		inline void output_segmentation_status_2(std::size_t const j, std::size_t const sample_count) const;
	};
}

#endif
