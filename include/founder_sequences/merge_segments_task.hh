/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef FOUNDER_SEQUENCES_MERGE_SEGMENTS_TASK_HH
#define FOUNDER_SEQUENCES_MERGE_SEGMENTS_TASK_HH

#include <atomic>
#include <founder_sequences/founder_sequences.hh>
#include <founder_sequences/segment_text.hh>
#include <founder_sequences/substring_copy_number.hh>
#include <lemon/full_graph.h>
#include <lemon/matching.h>
#include <vector>


namespace founder_sequences {
	
	class merge_segments_task;
	
	
	struct merge_segments_task_delegate
	{
		virtual ~merge_segments_task_delegate() {}
		virtual void task_did_finish(merge_segments_task &task) = 0;
	};
	
	
	class merge_segments_task final : public task
	{
	protected:
		typedef int32_t																weight_type;
		typedef lemon::FullBpGraph													graph_type;
		typedef int32_t																edge_cost_type;
		typedef std::vector <edge_cost_type>										edge_cost_vector_type;
		typedef graph_type::EdgeMap <edge_cost_type>								edge_cost_map_type;
		typedef lemon::MaxWeightedPerfectMatching <graph_type, edge_cost_map_type>	matching_type;	// No bipartite algorithms in Lemon 1.3.1.
		
	protected:
		merge_segments_task_delegate	*m_delegate{};
		segment_text_vector				*m_lhs{};
		segment_text_vector				*m_rhs{};
		matching_vector					*m_matching{};
		std::size_t						m_task_idx{};
		bipartite_set_scoring			m_bipartite_set_scoring_method{};
		
	public:
		merge_segments_task() = default;
		
		merge_segments_task(
			std::size_t task_idx,
			merge_segments_task_delegate &delegate,
			segment_text_vector &lhs,
			segment_text_vector &rhs,
			matching_vector &matching,
			bipartite_set_scoring bipartite_set_scoring_method
		):
			m_delegate(&delegate),
			m_lhs(&lhs),
			m_rhs(&rhs),
			m_matching(&matching),
			m_task_idx(task_idx),
			m_bipartite_set_scoring_method(bipartite_set_scoring_method)
		{
		}
		
		std::size_t task_index() const { return m_task_idx; }
		void execute() override;
		
	protected:
		inline std::size_t edge_id(
			graph_type const &graph,
			std::size_t const li,
			std::size_t const ri
		);
		
		void copy_edge_weight(
			graph_type const &graph,
			edge_cost_vector_type &edge_costs,
			std::size_t const li,
			std::size_t const ri,
			std::size_t const source_li,
			std::size_t const source_ri
		);
		
		void calculate_edge_weight(
			graph_type const &graph,
			edge_cost_vector_type &edge_costs,
			std::size_t const li,
			std::size_t const ri
		);
		
		weight_type find_minimum_cost_matching(
			graph_type const &graph,
			edge_cost_map_type const &edge_costs
		);
	};
}

#endif
