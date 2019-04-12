/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef FOUNDER_SEQUENCES_BIPARTITE_MATCHER_HH
#define FOUNDER_SEQUENCES_BIPARTITE_MATCHER_HH

#include <atomic>
#include <founder_sequences/founder_sequences.hh>
#include <founder_sequences/matcher.hh>
#include <founder_sequences/merge_segments_task.hh>
#include <founder_sequences/segment_text.hh>
#include <founder_sequences/substring_copy_number.hh>
#include <libbio/dispatch.hh>


namespace founder_sequences
{
	class bipartite_matcher;
	
	
	struct bipartite_matcher_delegate : public virtual matcher_delegate
	{
		virtual void matcher_did_finish(bipartite_matcher &matcher) = 0;
	};
	
	
	class bipartite_matcher final : public matcher, public merge_segments_task_delegate
	{
		// No point in substituting mutexes with atomic counters if they are not lock-free.
		static_assert(std::atomic_uint8_t::is_always_lock_free);
		
	protected:
		libbio::dispatch_ptr <dispatch_queue_t>		m_producer_queue;
		libbio::dispatch_ptr <dispatch_queue_t>		m_consumer_queue;
		std::vector <std::unique_ptr <task>>		m_tasks;
		std::vector <matching_vector> 				m_matchings;
		segment_text_matrix							m_segment_texts;
		
		// For merge_segments_tasks.
		std::vector <std::uint32_t>					m_segment_text_permutation;
		
		substring_copy_number_matrix const			*m_substrings_to_output{};
		bipartite_matcher_delegate					*m_delegate{};
		
	public:
		bipartite_matcher(
			bipartite_matcher_delegate &delegate,
			substring_copy_number_matrix const &substrings_to_output
		):
			m_producer_queue(delegate.producer_queue()),				// Copy.
			m_consumer_queue(delegate.consumer_queue()),				// Copy.
			m_substrings_to_output(&substrings_to_output),
			m_delegate(&delegate)
		{
		}
		
		void match() override;
		void output_segments(std::ostream &stream, sequence_vector const &sequences) override;
		void task_did_finish(merge_segments_task &task) override {}; // No-op.
		
	protected:
		void create_initial_permutation();
		void start_matching_tasks();
		void create_permutations_and_notify();
	};
}

#endif
