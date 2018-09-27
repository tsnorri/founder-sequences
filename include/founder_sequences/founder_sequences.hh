/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef FOUNDER_SEQUENCES_FOUNDER_SEQUENCES_HH
#define FOUNDER_SEQUENCES_FOUNDER_SEQUENCES_HH

#include <libbio/cxxcompat.hh>
#include <libbio/consecutive_alphabet.hh>
#include <libbio/dispatch_fn.hh>
#include <libbio/pbwt_context.hh>
#include <sdsl/int_vector.hpp>
#include <sdsl/rmq_support.hpp>
#include <vector>


namespace founder_sequences {
	
	enum class segment_joining : std::uint8_t {
		GREEDY = 0,
		BIPARTITE_MATCHING,
		RANDOM,
		PBWT_ORDER
	};
	
	enum class bipartite_set_scoring : std::uint8_t {
		SYMMETRIC_DIFFERENCE = 0,
		INTERSECTION
	};
	
	typedef std::span <std::uint8_t const>					sequence;
	typedef std::vector <std::uint32_t>						matching_vector;
	typedef std::vector <sequence>							sequence_vector;
	typedef libbio::consecutive_alphabet_as <std::uint8_t>	alphabet_type;
	
	typedef sdsl::int_vector <0>							permutation_vector;
	typedef std::vector <permutation_vector>				permutation_matrix;
	
	
	template <typename t_element>
	using vector_tpl = std::vector <t_element>;
	
	
	typedef libbio::pbwt::pbwt_context <
		sequence_vector,					/* sequence_vector */
		alphabet_type,						/* alphabet_type */
		sdsl::range_maximum_sct <>::type,	/* ci_rmq */
		sdsl::int_vector <32>,				/* string_index_vector */
		sdsl::int_vector <32>,				/* character_index_vector */
		sdsl::int_vector <32>,				/* count_vector */
		std::uint32_t						/* divergence_count_type */
	> pbwt_context;
	
	typedef libbio::pbwt::buffering_pbwt_context <
		sequence_vector,					/* sequence_vector */
		alphabet_type,						/* alphabet_type */
		sdsl::range_maximum_sct <>::type,	/* ci_rmq */
		std::uint32_t,						/* string_index */
		std::uint32_t,						/* character_index */
		std::uint32_t,						/* character_count */
		std::uint32_t,						/* divergence_count */
		vector_tpl							/* vector_tpl */
	> buffering_pbwt_context;
	
	
	typedef buffering_pbwt_context::sample_type		pbwt_sample_type;
	
	
	struct segmentation_context
	{
		virtual ~segmentation_context() {}
		virtual std::uint32_t max_segment_size() const = 0;
	};
	
	
	struct segmentation_context_delegate
	{
		virtual ~segmentation_context_delegate() {}
		virtual alphabet_type const &alphabet() const = 0;
		virtual sequence_vector const &sequences() const = 0;
		virtual bipartite_set_scoring bipartite_set_scoring_method() const = 0;
		virtual bool should_run_single_threaded() const = 0;
		
		virtual std::uint32_t sequence_count() const = 0;
		
		virtual std::ostream &sequence_output_stream() = 0;
		virtual std::ostream &segments_output_stream() = 0;
	};
	
	
	struct task
	{
		virtual ~task() {}
		virtual void execute() = 0;
	};
}

#endif
