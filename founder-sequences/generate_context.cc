/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <algorithm>
#include <boost/range/combine.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/range/combine.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/vector.hpp>
#include <experimental/iterator>
#include <founder_sequences/founder_sequences.hh>
#include <founder_sequences/generate_context.hh>
#include <founder_sequences/rmq.hh>
#include <libbio/algorithm.hh>
#include <libbio/consecutive_alphabet.hh>
#include <libbio/counting_sort.hh>
#include <libbio/fasta_reader.hh>
#include <libbio/line_reader.hh>
#include <libbio/sdsl_boost_serialization.hh>
#include <libbio/vector_source.hh>
#include <memory>
#include <vector>

namespace lb	= libbio;
namespace lsr	= libbio::sequence_reader;


namespace founder_sequences { namespace detail {
	
	std::size_t progress_indicator_gc_data_source::progress_step_max() const
	{
		return m_ctx->m_step_max;
	}
	
	
	std::size_t progress_indicator_gc_data_source::progress_current_step() const
	{
		return m_ctx->m_current_step.load(std::memory_order_relaxed);
	}
	
	
	void progress_indicator_lp_generate_traceback_data_source::progress_log_extra() const
	{
		auto const sample_count(m_context->current_pbwt_sample_count());
		
		std::cerr << ", " << sample_count;
		if (1 == sample_count)
			std::cerr << " sample";
		else
			std::cerr << " samples";
	}
}}


namespace founder_sequences {
	
	void generate_context::load_input(char const *input_path, lsr::input_format const input_file_format)
	{
		lb::log_time(std::cerr);
		std::cerr << "Loading the input…" << std::flush;
		
		lsr::read_input(input_path, input_file_format, m_sequence_container);
		m_sequence_container->to_spans(m_sequences);
		
		if (0 == m_sequences.size())
		{
			std::cerr << "\nThe input file contained no sequences." << std::endl;
			exit(EXIT_SUCCESS);
		}
		
		auto const seq_length(m_sequences.front().size());
		std::cerr << " length: " << seq_length << std::endl;
	}
	
	
	void generate_context::check_input() const
	{
		lb::log_time(std::cerr);
		std::cerr << "Checking the input…" << std::endl;

		// Check that all the vectors have equal lengths.
		auto const sequence_length(m_sequences.front().size());
		bool stop(false);
		std::size_t i(1);
		for (auto const &vec : m_sequences | ranges::view::drop(1))
		{
			if (vec.size() != sequence_length)
			{
				stop = true;
				std::cerr
					<< "The length of the sequence at index " << i << " was " << vec.size()
					<< " while that of the first one was " << sequence_length << '.' << std::endl;
			}
			++i;
		}
	
		if (stop)
			exit(EXIT_FAILURE);
	}
	
	
	void generate_context::generate_alphabet_and_continue()
	{
		auto const sequence_length(m_sequences.front().size());
		
		if (0 == m_pbwt_sample_rate)
		{
			libbio::log_time(std::cerr);
			std::cerr << "PBWT sampling shall not be done." << std::endl;
			m_pbwt_sample_rate = 1 + sequence_length;
		}
		else
		{
			auto const multiplier(m_pbwt_sample_rate);
			m_pbwt_sample_rate = std::ceil(multiplier * std::sqrt(sequence_length));
			libbio::log_time(std::cerr);
			std::cerr << "Using " << multiplier << "√n = " << m_pbwt_sample_rate << " as the sample rate." << std::endl;
		}

		libbio::log_time(std::cerr);
		std::cerr << "Generating a compressed alphabet…" << std::endl;

		m_current_step = 0;
		m_step_max = m_sequences.size();
		m_progress_indicator_data_source.reset(new detail::progress_indicator_gc_data_source(*this));
		
		dispatch_async(*m_parallel_queue, ^{
			lb::consecutive_alphabet_as_builder <std::uint8_t> builder;
			
			builder.init();
			std::size_t i(0);
			for (auto const &vec : m_sequences)
			{
				builder.prepare(vec);
				m_current_step.fetch_add(1, std::memory_order_relaxed);
			}
			builder.compress();
			
			using std::swap;
			swap(m_alphabet, builder.alphabet());

			dispatch_async(dispatch_get_main_queue(), ^{
				m_progress_indicator.end_logging_mt();
				m_current_step = 0;
				m_step_max = 0;
				calculate_segmentation(0, sequence_length);
			});
		});

		m_progress_indicator.log_with_progress_bar("\t", *m_progress_indicator_data_source);
	}
	
	
	void generate_context::calculate_segmentation_short_path(std::size_t const lb, std::size_t const rb)
	{
		segmentation_sp_context ctx(*this, lb, rb);
		ctx.process();
		std::cerr << "Outputting…" << std::endl;
		ctx.output_founders();
		
		if (m_segments_ostream_ptr)
			ctx.output_segments();
		
		// Finish.
		finish();
		lb::log_time(std::cerr);
		std::cerr << "Done." << std::endl;
		exit(EXIT_SUCCESS);
	}
	
	
	void generate_context::calculate_segmentation_long_path(std::size_t const lb, std::size_t const rb)
	{
		assert(m_segment_length < rb - lb);
		assert(0 == lb); // FIXME: handle ranges that don't start from zero.
		
		auto *ctx(new segmentation_lp_context(*this, m_parallel_queue, m_serial_queue)); // Uses callbacks, deleted in the final one.
		m_progress_indicator_data_source.reset(new detail::progress_indicator_lp_generate_traceback_data_source(*ctx));
		
		ctx->generate_traceback(lb, rb);
		m_progress_indicator.log_with_progress_bar("\t", *m_progress_indicator_data_source);
	}
	
	
	void generate_context::check_traceback_size(segmentation_context &ctx)
	{
		if (! (ctx.max_segment_size() < sequence_count()))
		{
			finish();
			std::cerr << "Unable to reduce the number of sequences; the maximum segment size is equal to the number of input sequences." << std::endl;
			exit(EXIT_FAILURE);
		}
	}
	
	
	void generate_context::context_did_finish_traceback(segmentation_sp_context &ctx)
	{
		check_traceback_size(ctx);
	}
	
	
	void generate_context::context_will_follow_traceback(segmentation_lp_context &ctx)
	{
		// Not main queue.
		
		// Update the progress indicator.
		m_progress_indicator.end_logging(); // Uses dispatch_sync.
		
		dispatch_async(dispatch_get_main_queue(), ^{
			lb::log_time(std::cerr);
			std::cerr << "Following the traceback…" << std::flush;
		});
		
	}
	
	
	void generate_context::context_did_finish_traceback(segmentation_lp_context &ctx, std::size_t const segment_count, std::size_t const max_segment_size)
	{
		assert(dispatch_get_current_queue() == dispatch_get_main_queue());
		
		std::cerr << " there were " << segment_count << " segments the maximum size of which was " << max_segment_size << '.' << std::endl;
		check_traceback_size(ctx);
		
		lb::log_time(std::cerr);
		std::cerr << "Updating the PBWT samples to traceback positions…" << std::endl;
		
		ctx.update_samples_to_traceback_positions();
	}
	
	
	void generate_context::context_will_start_update_samples_tasks(segmentation_lp_context &ctx)
	{
		// Not main queue.
		
		m_progress_indicator_data_source.reset(new detail::progress_indicator_lp_generic_data_source(ctx));
		m_progress_indicator.log_with_progress_bar("\t", *m_progress_indicator_data_source);
	}
	
	
	void generate_context::context_did_start_update_samples_tasks(segmentation_lp_context &ctx)
	{
		// Not main queue.
	}
	
	
	void generate_context::context_did_update_pbwt_samples_to_traceback_positions(segmentation_lp_context &ctx)
	{
		assert(dispatch_get_current_queue() == dispatch_get_main_queue());

		m_progress_indicator.end_logging_mt();
		
		lb::log_time(std::cerr);
		std::cerr << "Reducing the number of segments…" << std::endl;
		
		ctx.find_segments_greedy();
	}
	
	
	void generate_context::context_will_merge_segments(segmentation_lp_context &ctx)
	{
		// Not main queue.
		
		m_progress_indicator_data_source.reset(new detail::progress_indicator_lp_generic_data_source(ctx));
		m_progress_indicator.log_with_progress_bar("\t", *m_progress_indicator_data_source);
	}
	
	
	void generate_context::context_did_merge_segments(segmentation_lp_context &ctx, segmentation_container &&container)
	{
		assert(dispatch_get_current_queue() == dispatch_get_main_queue());
		
		m_progress_indicator.end_logging_mt();
		
		// Context no longer needed, deallocate.
		ctx.cleanup();
		
		if (m_segmentation_ostream.is_open())
			save_segmentation_to_file(container);
		
		if (running_mode::GENERATE_FOUNDERS == m_running_mode)
			join_segments_and_output(std::move(container));
		else
			finish_lp();
	}
	
	
	void generate_context::load_segmentation_from_file(segmentation_container &container)
	{
		assert(dispatch_get_current_queue() == dispatch_get_main_queue());
		assert(m_segmentation_istream.is_open());
		
		lb::log_time(std::cerr);
		std::cerr << "Loading the segmentation…" << std::flush;
		std::string stored_input_path;
		boost::archive::text_iarchive archive(m_segmentation_istream);
		archive >> stored_input_path;
		archive >> m_segment_length;
		archive >> m_alphabet;
		archive >> container;
		
		std::cerr << " done.\n";
		std::cerr << "\tStored input path: '" << stored_input_path << "'\n";
		std::cerr << "\tSegment length bound: " << m_segment_length << std::endl;
	}

	
	void generate_context::save_segmentation_to_file(segmentation_container const &container)
	{
		assert(dispatch_get_current_queue() == dispatch_get_main_queue());
		assert(m_segmentation_ostream.is_open());
		
		lb::log_time(std::cerr);
		std::cerr << "Saving the segmentation…" << std::endl;
		
		boost::archive::text_oarchive archive(m_segmentation_ostream);
		archive << m_sequence_container->path();
		archive << m_segment_length;
		archive << m_alphabet;
		archive << container;
	}
	
	
	void generate_context::join_segments_and_output(segmentation_container &&container)
	{
		assert(dispatch_get_current_queue() == dispatch_get_main_queue());
		
		lb::log_time(std::cerr);
		std::cerr << "Joining the remaining segments…" << std::endl;
		
		auto *ctx(new join_context(*this, m_parallel_queue, m_serial_queue, std::move(container), m_random_seed)); // Uses callbacks, deleted in the final one.
		ctx->join_segments_and_output(m_segment_joining_method);
	}
	
	
	void generate_context::context_will_output_founders(join_context &ctx)
	{
		assert(dispatch_get_current_queue() == dispatch_get_main_queue());
		
		lb::log_time(std::cerr);
		std::cerr << "Outputting the founders…" << std::endl;
	}
	
	
	void generate_context::context_did_output_founders(join_context &ctx)
	{
		assert(dispatch_get_current_queue() == dispatch_get_main_queue());
		
		if (m_segments_ostream_ptr)
		{
			lb::log_time(std::cerr);
			std::cerr << "Outputting the segments…" << std::endl;
			ctx.output_segments(m_segment_joining_method);
		}
		
		ctx.cleanup();
		finish_lp();
	}
	
	
	void generate_context::finish_lp()
	{
		assert(dispatch_get_current_queue() == dispatch_get_main_queue());
		
		// Finish.
		finish();
		lb::log_time(std::cerr);
		std::cerr << "Done." << std::endl;
		exit(EXIT_SUCCESS);
	}
	
	
	void generate_context::calculate_segmentation(std::size_t const lb, std::size_t const rb)
	{
		assert(dispatch_get_current_queue() == dispatch_get_main_queue());
		
		lb::log_time(std::cerr);
		std::cerr << "Calculating the segmentation…" << std::endl;
		
		if (rb - lb < 2 * m_segment_length)
			calculate_segmentation_short_path(lb, rb);
		else
			calculate_segmentation_long_path(lb, rb);
	}
	
	
	void generate_context::prepare(
		char const *segmentation_input_path,
		char const *segmentation_output_path,
		char const *output_founders_path,
		char const *output_segments_path
	)
	{
		if (m_use_single_thread)
		{
			lb::dispatch_ptr <dispatch_queue_t> queue(dispatch_queue_create("fi.iki.tsnorri.worker-queue", DISPATCH_QUEUE_SERIAL), false);
			m_parallel_queue = queue;
			m_serial_queue = queue;
		}
		else
		{
			m_parallel_queue.reset(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), true);
			m_serial_queue.reset(dispatch_queue_create("fi.iki.tsnorri.processing-queue", DISPATCH_QUEUE_SERIAL), false);
		}
		
		if (segmentation_input_path)
			lb::open_file_for_reading(segmentation_input_path, m_segmentation_istream);
		else if (segmentation_output_path)
			lb::open_file_for_writing(segmentation_output_path, m_segmentation_ostream, lb::writing_open_mode::CREATE);
		
		if (output_founders_path)
			lb::open_file_for_writing(output_founders_path, m_founders_ostream, lb::writing_open_mode::CREATE);
		
		if (output_segments_path)
		{
			if ('-' == output_segments_path[0] && '\0' == output_segments_path[1])
				m_segments_ostream_ptr = &std::cout;
			else
			{
				lb::open_file_for_writing(output_segments_path, m_segments_ostream, lb::writing_open_mode::CREATE);
				m_segments_ostream_ptr = &m_segments_ostream;
			}
		}

		if (m_progress_indicator.is_stderr_interactive())
			m_progress_indicator.install();
	}
	
	
	void generate_context::load_and_generate(
		char const *input_path,
		lsr::input_format const input_file_format
	)
	{
		load_input(input_path, input_file_format);
		check_input();
		
		if (m_segmentation_istream.is_open())
		{
			segmentation_container container;
			load_segmentation_from_file(container);
			join_segments_and_output(std::move(container));
		}
		else
		{
			generate_alphabet_and_continue();
		}
	}
	
	
	void generate_context::finish()
	{
		m_progress_indicator.uninstall();
		cleanup();
	}
}

