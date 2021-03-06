/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef FOUNDER_SEQUENCES_GENERATE_CONTEXT_HH
#define FOUNDER_SEQUENCES_GENERATE_CONTEXT_HH

#include <founder_sequences/founder_sequences.hh>
#include <founder_sequences/join_context.hh>
#include <founder_sequences/segmentation_dp_arg.hh>
#include <founder_sequences/segmentation_lp_context.hh>
#include <founder_sequences/segmentation_sp_context.hh>
#include <libbio/dispatch.hh>
#include <libbio/file_handling.hh>
#include <libbio/sequence_reader/sequence_reader.hh>


namespace founder_sequences {
	class generate_context;
}


namespace founder_sequences { namespace detail {

	struct progress_indicator_data_source : public libbio::progress_indicator_delegate {};

	class progress_indicator_gc_data_source : public progress_indicator_data_source
	{
	protected:
		generate_context const	*m_ctx{};

	public:
		progress_indicator_gc_data_source() = default;
		progress_indicator_gc_data_source(generate_context const &ctx): m_ctx(&ctx) {}
		std::size_t progress_step_max() const override;
		std::size_t progress_current_step() const override;
		void progress_log_extra() const override {}
	};
	
	class progress_indicator_lp_data_source : public progress_indicator_data_source
	{
	protected:
		segmentation_lp_context	*m_context{};
		
	public:
		progress_indicator_lp_data_source() = default;
		progress_indicator_lp_data_source(segmentation_lp_context &ctx):
			m_context(&ctx)
		{
		}
		
		std::size_t progress_step_max() const override { return m_context->step_max(); }
		std::size_t progress_current_step() const override { return m_context->current_step(); }
	};
	
	struct progress_indicator_lp_generate_traceback_data_source final : public progress_indicator_lp_data_source
	{
		using progress_indicator_lp_data_source::progress_indicator_lp_data_source;
		void progress_log_extra() const override;
	};
	
	struct progress_indicator_lp_generic_data_source final : public progress_indicator_lp_data_source
	{
		using progress_indicator_lp_data_source::progress_indicator_lp_data_source;
		void progress_log_extra() const override {}
	};
}}


namespace founder_sequences {
	
	class generate_context final : public segmentation_lp_context_delegate, public segmentation_sp_context_delegate, public join_context_delegate
	{
		friend class detail::progress_indicator_gc_data_source;
		
	protected:
		libbio::dispatch_ptr <dispatch_queue_t>							m_parallel_queue;
		libbio::dispatch_ptr <dispatch_queue_t>							m_serial_queue;
		
		std::unique_ptr <libbio::sequence_reader::sequence_container>	m_sequence_container;
		sequence_vector													m_sequences;
		alphabet_type													m_alphabet;
		
		libbio::file_istream											m_segmentation_istream;
		libbio::file_ostream											m_segmentation_ostream;
		libbio::file_ostream											m_founders_ostream;
		libbio::file_ostream											m_segments_ostream;
		std::ostream													*m_segments_ostream_ptr{};
		
		libbio::progress_indicator										m_progress_indicator;
		std::unique_ptr <detail::progress_indicator_data_source>		m_progress_indicator_data_source;

		std::atomic_uint32_t											m_current_step{};
		std::atomic_uint32_t											m_step_max{};
		
		std::size_t														m_segment_length{};
		std::uint64_t													m_pbwt_sample_rate{};
		std::uint_fast32_t												m_random_seed{};
		running_mode													m_running_mode{};
		segment_joining													m_segment_joining_method{};
		bipartite_set_scoring											m_bipartite_set_scoring{};
		bool															m_use_single_thread{false};
	
	public:
		generate_context(
			running_mode const mode,
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
			m_running_mode(mode),
			m_segment_joining_method(segment_joining_method),
			m_bipartite_set_scoring(bipartite_set_scoring),
			m_use_single_thread(use_single_thread)
		{
		}
		
		generate_context(generate_context const &) = delete;
		generate_context(generate_context &&) = delete;
		
		sequence_vector const &sequences() const override { return m_sequences; }
		std::uint32_t sequence_count() const override { return m_sequences.size(); }
		alphabet_type const &alphabet() const override { return m_alphabet; }
		std::size_t segment_length() const override { return m_segment_length; }
		std::uint64_t pbwt_sample_rate() const override { return m_pbwt_sample_rate; }
		std::ostream &sequence_output_stream() override { return (m_founders_ostream.is_open() ? m_founders_ostream : std::cout); }
		std::ostream &segments_output_stream() override { return *m_segments_ostream_ptr; }
		bipartite_set_scoring bipartite_set_scoring_method() const override { return m_bipartite_set_scoring; }
		bool should_run_single_threaded() const override { return m_use_single_thread; }

		void context_did_finish_traceback(segmentation_sp_context &ctx) override;
		
		void context_will_follow_traceback(segmentation_lp_context &ctx) override;
		void context_did_finish_traceback(segmentation_lp_context &ctx, std::size_t const segment_count, std::size_t const max_segment_size) override;
		void context_will_start_update_samples_tasks(segmentation_lp_context &ctx) override;
		void context_did_start_update_samples_tasks(segmentation_lp_context &ctx) override;
		void context_did_update_pbwt_samples_to_traceback_positions(segmentation_lp_context &ctx) override;
		void context_will_merge_segments(segmentation_lp_context &ctx) override;
		void context_did_merge_segments(segmentation_lp_context &ctx, segmentation_container &&container) override;
		
		void join_segments_and_output(segmentation_container &&container);
		
		void context_will_output_founders(join_context &ctx) override;
		void context_did_output_founders(join_context &ctx) override;
		
		void prepare(
			char const *segmentation_input_path,
			char const *segmentation_output_path,
			char const *output_founders_path,
			char const *output_segments_path
		);
		void load_and_generate(
			char const *input_path,
			libbio::sequence_reader::input_format const input_file_format
		);
		void finish();
		
		// For debugging.
		void print_segment_texts() const;
		
	protected:
		void load_input(char const *input_path, libbio::sequence_reader::input_format const input_file_format);
		void check_input() const;
		void generate_alphabet_and_continue();
		void generate_founders(std::size_t const lb, std::size_t const rb);
	
		void calculate_segmentation(std::size_t const lb, std::size_t const rb);
		void calculate_segmentation_short_path(std::size_t const lb, std::size_t const rb);
		void calculate_segmentation_long_path(std::size_t const lb, std::size_t const rb);
		
		void check_traceback_size(segmentation_context &ctx);
		
		void load_segmentation_from_file(segmentation_container &container);
		void save_segmentation_to_file(segmentation_container const &container);
		
		void finish_lp();
		void cleanup() { delete this; }
	};
}

#endif
