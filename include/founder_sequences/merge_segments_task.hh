/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef FOUNDER_SEQUENCES_MERGE_SEGMENTS_TASK_HH
#define FOUNDER_SEQUENCES_MERGE_SEGMENTS_TASK_HH

#include <dispatch/dispatch.h>
#include <founder_sequences/founder_sequences.hh>
#include <lemon/full_graph.h>
#include <lemon/matching.h>
#include <vector>


namespace founder_sequences {
	
	class merge_segments_task final
	{
	protected:
		typedef int32_t																weight_type;
		typedef lemon::FullBpGraph													graph_type;
		typedef int32_t																edge_cost_type;
		typedef std::vector <edge_cost_type>										edge_cost_vector_type;
		typedef graph_type::EdgeMap <edge_cost_type>								edge_cost_map_type;
		typedef lemon::MaxWeightedPerfectMatching <graph_type, edge_cost_map_type>	matching_type;	// No bipartite algorithms in Lemon 1.3.1.
		
	protected:
		std::size_t			m_task_idx{0};
		segment_text_vector	*m_lhs{nullptr};
		segment_text_vector	*m_rhs{nullptr};
		
	public:
		merge_segments_task() = default;
		
		merge_segments_task(
			std::size_t task_idx,
			segment_text_vector &lhs,
			segment_text_vector &rhs
		):
			m_task_idx(task_idx),
			m_lhs(&lhs),
			m_rhs(&rhs)
		{
		}
		
		void execute(dispatch_queue_t queue, matching_vector &matchings);
		
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
			edge_cost_map_type const &edge_costs,
			matching_vector &matchings
		);
	};
}

#endif
