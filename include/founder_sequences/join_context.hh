/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef FOUNDER_SEQUENCES_JOIN_CONTEXT_HH
#define FOUNDER_SEQUENCES_JOIN_CONTEXT_HH

#include <founder_sequences/bipartite_matcher.hh>
#include <founder_sequences/greedy_matcher.hh>
#include <founder_sequences/segmentation_container.hh>
#include <founder_sequences/segmentation_context.hh>
#include <libbio/dispatch_fn.hh>


namespace founder_sequences {
	
	class join_context;
	
	
	struct join_context_delegate : public virtual segmentation_context_delegate
	{
		virtual void context_will_output_founders(join_context &ctx) = 0;
		virtual void context_did_output_founders(join_context &ctx) = 0;
	};
	
	
	class join_context final : public bipartite_matcher_delegate, public greedy_matcher_delegate
	{
	protected:
		libbio::dispatch_ptr <dispatch_queue_t>				m_producer_queue;	// May be parallel.
		libbio::dispatch_ptr <dispatch_queue_t>				m_consumer_queue;	// Needs to be serial.
		
		segmentation_container								m_segmentation_container{};
		std::uint_fast32_t									m_random_seed{};
		
		// For matching.
		substring_copy_number_matrix						m_substring_copy_numbers;
		permutation_matrix									m_permutations;
		std::unique_ptr <matcher>							m_matcher;
		std::uint32_t										m_permutation_max{};
		std::uint8_t										m_permutation_bits_needed{};
		
		join_context_delegate								*m_delegate{};
		
	public:
		join_context(
			join_context_delegate &delegate,
			libbio::dispatch_ptr <dispatch_queue_t> &producer_queue,
			libbio::dispatch_ptr <dispatch_queue_t> &consumer_queue,
			segmentation_container &&segmentation_container,
			std::uint_fast32_t const random_seed
		):
			m_producer_queue(producer_queue),
			m_consumer_queue(consumer_queue),
			m_segmentation_container(std::move(segmentation_container)),
			m_random_seed(random_seed),
			m_delegate(&delegate)
		{
		}
		
		void cleanup() { delete this; }
		
		void join_segments_and_output(segment_joining const seg_joining);
		void output_segments(segment_joining const seg_joining) const;	// FIXME: segment_joining not actually needed here if random and PBWT order output are moved to matcher classes.
		
		libbio::dispatch_ptr <dispatch_queue_t> producer_queue() const override { return m_producer_queue; }
		libbio::dispatch_ptr <dispatch_queue_t> consumer_queue() const override { return m_consumer_queue; }
		std::uint32_t sequence_count() const override { return m_delegate->sequence_count(); }
		std::uint32_t max_segment_size() const override { return m_segmentation_container.max_segment_size; }
		std::vector <pbwt_sample_type> const &pbwt_samples() const override { return m_segmentation_container.reduced_pbwt_samples; }
		segmentation_traceback_vector const &reduced_traceback() const override { return m_segmentation_container.reduced_traceback; }
		permutation_matrix &permutations() override { return m_permutations; }
		bipartite_set_scoring bipartite_set_scoring_method() const override { return m_delegate->bipartite_set_scoring_method(); }
		
		void matcher_did_finish(bipartite_matcher &matcher) override;
		void matcher_did_finish(greedy_matcher &matcher) override;
		
	protected:
		void make_cumulative_sum(substring_copy_number_vector &vec) const;
		
		void init_permutations();
		
		void join_greedy();
		void join_with_bipartite_matching();
		void join_random_order_and_output();
		void join_pbwt_order_and_output() const;
		
		void output_in_permutation_order() const;
		
		void output_substring(
			std::ostream &ostream,
			std::size_t const sequence_idx,
			std::size_t const text_pos,
			std::size_t const text_length,
			sequence_vector const &sequences
		) const;
		
		void output_gaps(
			std::ostream &ostream,
			std::size_t const text_pos,
			std::size_t const text_length
		) const;
	};
}

#endif
