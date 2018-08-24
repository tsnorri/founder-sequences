/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <dispatch/dispatch.h>
#include <founder_sequences/merge_segments_task.hh>
#include <libbio/algorithm.hh>
#include <libbio/assert.hh>


namespace lb	= libbio;


namespace founder_sequences {
	
	std::size_t merge_segments_task::edge_id(
		graph_type const &graph,
		std::size_t const li,
		std::size_t const ri
	)
	{
		auto const lhs(graph.redNode(li));
		auto const rhs(graph.blueNode(ri));
		auto const edge(graph.edge(lhs, rhs));
		auto const edge_id(graph_type::id(edge));
		return edge_id;
	}
	
	
	void merge_segments_task::calculate_edge_weight(
		graph_type const &graph,
		edge_cost_vector_type &edge_costs,
		std::size_t const li,
		std::size_t const ri
	)
	{
		auto const &lhs((*m_lhs)[li]);
		auto const &rhs((*m_rhs)[ri]);
		
		auto const opposite_weight(lb::set_symmetric_difference_size(
			lhs.sequence_indices.cbegin(),
			lhs.sequence_indices.cend(),
			rhs.sequence_indices.cbegin(),
			rhs.sequence_indices.cend()
		));

		// Try to make sure that the weight is within data type limits.
		weight_type const weight(-opposite_weight);
		libbio_always_assert(-weight == opposite_weight);
		
		// Store the weight.
		auto const lhs_node(graph.redNode(li));
		auto const rhs_node(graph.blueNode(ri));
		auto const edge(graph.edge(lhs_node, rhs_node));
		auto const edge_id(graph_type::id(edge));
		libbio_always_assert(edge_id < edge_costs.size());
		edge_costs[edge_id] = weight;
	}
	
	
	void merge_segments_task::copy_edge_weight(
		graph_type const &graph,
		edge_cost_vector_type &edge_costs,
		std::size_t const li,
		std::size_t const ri,
		std::size_t const source_li,
		std::size_t const source_ri
	)
	{
		// Fetch and store the weight.
		auto const source_eid(edge_id(graph, source_li, source_ri));
		auto const eid(edge_id(graph, li, ri));
		libbio_always_assert(eid < edge_costs.size());
		libbio_always_assert(source_eid < edge_costs.size());
		edge_costs[eid] = edge_costs[source_eid];
	}
	
	
	auto merge_segments_task::find_minimum_cost_matching(
		graph_type const &graph,
		edge_cost_map_type const &edge_costs,
		matching_vector &matchings
	) -> weight_type
	{
		matching_type matching(graph, edge_costs);
		matching.run();

		// Iterate red nodes in the graph, find their mates and store the pairs in target_paths.
		auto const matching_weight(matching.matchingWeight());
		for (std::size_t i(0), count(graph.redNum()); i < count; ++i)
		{
			auto const lhs(graph.redNode(i));
			libbio_always_assert(lhs != lemon::INVALID);
			auto const node(matching.mate(lhs));
			libbio_always_assert(node != lemon::INVALID);
			auto const rhs(graph.asBlueNode(node));
			libbio_always_assert(rhs != lemon::INVALID);
			
			auto const li(graph.index(lhs));
			auto const ri(graph.index(rhs));
			libbio_always_assert(li < matchings.size());
			libbio_always_assert(ri < matchings.size());
			matchings[li] = ri;
		}
		
		return matching_weight;
	}
	
	
	void merge_segments_task::execute(matching_vector &matchings)
	{
		auto const path_count(m_lhs->size());
		auto const rhs_path_count(m_rhs->size());
		libbio_always_assert(rhs_path_count == path_count);
		//dispatch_ptr <dispatch_group_t> group(dispatch_group_create());
		
		// Use a std::vector since contents of its contained objects in different elements
		// may be modified concurrently without data races.
		graph_type graph(path_count, path_count);
		edge_cost_map_type edge_cost_map(graph);
		
		{
			std::size_t const edge_costs_size(1 + graph.maxEdgeId());
			edge_cost_vector_type edge_costs(edge_costs_size, 0);
			for (std::size_t i(0); i < path_count; ++i)
			{
				for (std::size_t j(0); j < path_count; ++j)
				{
					auto const &lhs((*m_lhs)[i]);
					auto const &rhs((*m_rhs)[j]);
					if (! (lhs.is_copied() || rhs.is_copied()))
					{
						//dispatch_group_async(*group, queue, ^{
							calculate_edge_weight(graph, edge_costs, i, j);
						//});
					}
				}
			}
			
			//dispatch_group_wait(*group, DISPATCH_TIME_FOREVER);
			
			// Copy the remaining edge weights.
			for (std::size_t i(0); i < path_count; ++i)
			{
				for (std::size_t j(0); j < path_count; ++j)
				{
					auto const &lhs((*m_lhs)[i]);
					auto const &rhs((*m_rhs)[j]);
					if (lhs.is_copied() || rhs.is_copied())
						copy_edge_weight(graph, edge_costs, i, j, lhs.row_number(i), rhs.row_number(j));
				}
			}
			
			// Copy the values to edge_cost_map.
			std::size_t i(0);
			for (auto const weight : edge_costs)
			{
				auto const edge(graph_type::edgeFromId(i));
				edge_cost_map.set(edge, weight);
				++i;
			}
		}
		
		libbio_always_assert(0 == matchings.size());
		matchings.resize(
			path_count,
			std::numeric_limits <matching_vector::value_type>::max()
		);
		auto const matching_weight(find_minimum_cost_matching(graph, edge_cost_map, matchings));

		m_done = true;
	}
}
