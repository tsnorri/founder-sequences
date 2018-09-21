/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef FOUNDER_SEQUENCES_GENERATE_CONTEXT_HH
#define FOUNDER_SEQUENCES_GENERATE_CONTEXT_HH

#include <founder_sequences/founder_sequences.hh>
#include <founder_sequences/segmentation_dp_arg.hh>
#include <founder_sequences/segmentation_lp_context.hh>
#include <founder_sequences/segmentation_sp_context.hh>
#include <libbio/dispatch_fn.hh>
#include <libbio/file_handling.hh>
#include <libbio/sequence_reader/sequence_reader.hh>


namespace founder_sequences {
	
	class generate_context final : public segmentation_lp_context_delegate, public segmentation_sp_context_delegate
	{
	protected:
		libbio::dispatch_ptr <dispatch_queue_t>							m_parallel_queue;
		libbio::dispatch_ptr <dispatch_queue_t>							m_serial_queue;
		
		std::unique_ptr <libbio::sequence_reader::sequence_container>	m_sequence_container;
		sequence_vector													m_sequences;
		alphabet_type													m_alphabet;
		
		libbio::file_ostream											m_founders_ostream;
		libbio::file_ostream											m_segments_ostream;
		std::ostream													*m_segments_ostream_ptr{};
		
		std::size_t														m_segment_length{};
		std::uint64_t													m_pbwt_sample_rate{};
		std::uint_fast32_t												m_random_seed{};
		segment_joining													m_segment_joining_method{};
		bipartite_set_scoring											m_bipartite_set_scoring{};
		bool															m_use_single_thread{false};
	
	public:
		generate_context(
			std::size_t const segment_length,
			segment_joining const segment_joining_method,
			bipartite_set_scoring const bipartite_set_scoring,
			std::uint64_t const pbwt_sample_rate,
			std::uint_fast32_t const random_seed,
			bool const use_single_thread
		):
			m_segment_length(segment_length),
			m_pbwt_sample_rate(pbwt_sample_rate),
			m_random_seed(random_seed),
			m_segment_joining_method(segment_joining_method),
			m_bipartite_set_scoring(bipartite_set_scoring),
			m_use_single_thread(use_single_thread)
		{
		}
	
		generate_context(generate_context const &) = delete;
		generate_context(generate_context &&) = delete;
	
		sequence_vector const &sequences() const override { return m_sequences; }
		alphabet_type const &alphabet() const override { return m_alphabet; }
		std::size_t segment_length() const override { return m_segment_length; }
		std::uint64_t pbwt_sample_rate() const override { return m_pbwt_sample_rate; }
		std::ostream &sequence_output_stream() override { return m_founders_ostream; }
		std::ostream &segments_output_stream() override { return *m_segments_ostream_ptr; }
		bipartite_set_scoring bipartite_set_scoring_method() const override { return m_bipartite_set_scoring; }

		void context_did_finish_traceback(segmentation_sp_context &ctx) override;

		void context_did_finish_traceback(segmentation_lp_context &ctx) override;
		void context_did_update_pbwt_samples_to_traceback_positions(segmentation_lp_context &ctx) override;
		void context_did_output_founders(segmentation_lp_context &ctx) override;

		void prepare(char const *output_founders_path, char const *output_segments_path);
		void load_and_generate(
			char const *input_path,
			libbio::sequence_reader::input_format const input_file_format
		);
		void cleanup() { delete this; }
	
		// For debugging.
		void print_segment_texts() const;
	
	protected:
		void load_input(char const *input_path, libbio::sequence_reader::input_format const input_file_format);
		void check_input() const;
		void generate_alphabet();
		void generate_founders(std::size_t const lb, std::size_t const rb);
	
		template <typename t_builder>
		void generate_alphabet(t_builder &builder);
	
		void calculate_segmentation(std::size_t const lb, std::size_t const rb);
		void calculate_segmentation_short_path(std::size_t const lb, std::size_t const rb);
		void calculate_segmentation_long_path(std::size_t const lb, std::size_t const rb);
		
		void check_traceback_size(segmentation_context &ctx);
	};
	
	
	template <typename t_builder>
	void generate_context::generate_alphabet(t_builder &builder)
	{
		libbio::log_time(std::cerr);
		std::cerr << "Generating a compressed alphabet…" << std::endl;
		
		builder.init();
		for (auto const &vec : m_sequences)
			builder.prepare(vec);
		builder.compress();
		
		using std::swap;
		swap(m_alphabet, builder.alphabet());
	}
}

#endif
