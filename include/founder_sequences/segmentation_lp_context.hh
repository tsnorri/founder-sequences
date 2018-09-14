/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef FOUNDER_SEQUENCES_SEGMENTATION_LP_CONTEXT_HH
#define FOUNDER_SEQUENCES_SEGMENTATION_LP_CONTEXT_HH

#include <founder_sequences/founder_sequences.hh>
#include <founder_sequences/segmentation_dp_arg.hh>
#include <founder_sequences/substring_copy_number.hh>
#include <founder_sequences/update_pbwt_task.hh>
#include <libbio/dispatch_fn.hh>


namespace founder_sequences { namespace detail {

	struct substring_index_pair
	{
		std::uint32_t	lhs_idx{};
		std::uint32_t	rhs_idx{};
		std::uint32_t	count{};
		
		substring_index_pair() = default;
		
		substring_index_pair(std::uint32_t lhs_idx_, std::uint32_t rhs_idx_, std::uint32_t count_):
			lhs_idx(lhs_idx_),
			rhs_idx(rhs_idx_),
			count(count_)
		{
		}
		
		substring_index_pair(std::uint32_t lhs_idx_, std::uint32_t rhs_idx_):
			lhs_idx(lhs_idx_),
			rhs_idx(rhs_idx_)
		{
		}
	};
	
	inline bool operator==(substring_index_pair const &lhs, substring_index_pair const &rhs)
	{
		return lhs.lhs_idx == rhs.rhs_idx && lhs.rhs_idx == rhs.rhs_idx;
	}
	
	inline std::ostream &operator<<(std::ostream &os, substring_index_pair const &pair)
	{
		os << '(' << pair.lhs_idx << ", " << pair.rhs_idx << "; " << pair.count << ')';
		return os;
	}
}}


namespace founder_sequences {
	
	void calculate_segmentation_lp_dp_arg(
		typename buffering_pbwt_context::divergence_count_list const &divergence_value_counts,
		segmentation_traceback_vector const &segmentation_traceback_dp,
		segmentation_traceback_vector_rmq const &segmentation_traceback_rmq,
		std::size_t const seq_count,
		std::size_t const segment_length,
		std::size_t const lb,	// Inclusive.
		std::size_t const text_pos,
		segmentation_dp_arg &min_arg
	);
	
	
	class segmentation_lp_context;
	
	
	struct segmentation_lp_context_delegate : public virtual segmentation_context_delegate
	{
		virtual std::size_t segment_length() const = 0;
		virtual alphabet_type const &alphabet() const = 0;
		virtual sequence_vector const &sequences() const = 0;
		virtual void context_did_finish_traceback(segmentation_lp_context &ctx) = 0;
		virtual void context_did_update_pbwt_samples_to_traceback_positions(segmentation_lp_context &ctx) = 0;
		virtual void context_did_output_founders(segmentation_lp_context &ctx) = 0;
	};
	
	
	class segmentation_lp_context final : public segmentation_context
	{
	protected:
		typedef std::vector <std::size_t>					text_position_vector;
		typedef buffering_pbwt_context::sample_type			pbwt_sample;
		typedef std::vector <substring_copy_number>			substring_copy_number_vector;
		typedef std::vector <substring_copy_number_vector>	substring_copy_number_matrix;
		typedef std::vector <sdsl::int_vector <0>>			permutation_vector;

	protected:
		buffering_pbwt_context								m_pbwt_ctx;
		segmentation_traceback_vector						m_segmentation_traceback_res;
		segmentation_traceback_vector						m_segmentation_traceback_dp;
		segmentation_traceback_vector_rmq					m_segmentation_traceback_dp_rmq;
		std::uint32_t										m_max_segment_size{};
		std::uint_fast32_t									m_random_seed{};
		
		libbio::dispatch_ptr <dispatch_queue_t>				m_producer_queue;	// May be parallel.
		libbio::dispatch_ptr <dispatch_queue_t>				m_consumer_queue;	// Needs to be serial.
		
		// For updating the PBWT samples.
		std::vector <std::unique_ptr <update_pbwt_task>>	m_update_pbwt_tasks;
		libbio::dispatch_ptr <dispatch_group_t>				m_update_samples_group;
		
		// For finding the reduced segmentation points.
		std::vector <pbwt_sample>							m_reduced_pbwt_samples;
		segmentation_traceback_vector						m_reduced_traceback;
		
		segmentation_lp_context_delegate					*m_delegate{};
		
	public:
		segmentation_lp_context(
			segmentation_lp_context_delegate &delegate,
			libbio::dispatch_ptr <dispatch_queue_t> &producer_queue,
			libbio::dispatch_ptr <dispatch_queue_t> &consumer_queue,
			std::uint_fast32_t const random_seed
		):
			m_pbwt_ctx(8, delegate.sequences(), delegate.alphabet()),
			m_random_seed(random_seed),
			m_producer_queue(producer_queue),
			m_consumer_queue(consumer_queue),
			m_delegate(&delegate)
		{
		}
		
		segmentation_traceback_vector &segmentation_traceback() { return m_segmentation_traceback_res; }
		std::uint32_t max_segment_size() const override { return m_max_segment_size; }
		std::uint32_t sequence_count() const override { return m_pbwt_ctx.size(); }
		
		void cleanup() { delete this; }
		
		void generate_traceback(std::size_t lb, std::size_t rb);
		void update_samples_to_traceback_positions();
		void find_segments_greedy();
		void join_segments_and_output(std::ostream &os, sequence_vector const &sequences, segment_joining const seg_joining);
		
	protected:
		void generate_traceback_part_2(std::size_t const lb, std::size_t const rb);
		void generate_traceback_part_3(std::size_t const lb, std::size_t const rb);
		void generate_traceback_part_4(std::size_t const lb, std::size_t const rb);
		void follow_traceback();
		
		void start_update_sample_task(
			std::size_t const lb,
			pbwt_sample &&sample,
			text_position_vector &&right_bounds
		);
		
		void make_cumulative_sum(substring_copy_number_vector &vec) const;
		
		std::pair <std::uint32_t, std::uint8_t> init_permutations(permutation_vector &permutations) const;
		
		void join_greedy_and_output(
			std::ostream &stream,
			substring_copy_number_matrix const &substrings_to_output,
			sequence_vector const &sequences
		) const;
		
		void join_with_bipartite_matching_and_output(
			std::ostream &stream,
			substring_copy_number_matrix const &substrings_to_output,
			sequence_vector const &sequences
		);
		
		void join_random_order_and_output(
			std::ostream &stream,
			substring_copy_number_matrix const &substrings_to_output,
			sequence_vector const &sequences
		) const;
		
		void join_pbwt_order_and_output(
			std::ostream &stream,
			substring_copy_number_matrix const &substrings_to_output,
			sequence_vector const &sequences
		) const;
		
		void output_in_permutation_order(
			std::ostream &stream,
			permutation_vector const &permutations,
			sequence_vector const &sequences,
			std::uint32_t const permutation_max
		) const;
		
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
		
		std::pair <std::size_t, std::size_t>
		greedy_create_index_pairs(
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
		) const;
		
		void greedy_create_matching(
			std::vector <detail::substring_index_pair> const &index_pairs,
			std::uint64_t const matching_max,
			sdsl::int_vector <0> &lhs_matching,
			sdsl::int_vector <0> &rhs_matching
		) const;


		inline void output_segmentation_status(std::size_t const j) const;
		inline void output_segmentation_status_mq(std::size_t const j) const;
	};
}

#endif
