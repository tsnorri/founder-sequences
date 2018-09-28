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
#include <libbio/vector_source.hh>
#include <memory>
#include <vector>

namespace lb	= libbio;
namespace lsr	= libbio::sequence_reader;


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
	
	
	void generate_context::generate_alphabet()
	{
		if (m_use_single_thread)
		{
			lb::consecutive_alphabet_as_builder <std::uint8_t> builder;
			generate_alphabet(builder);
		}
		else
		{
			lb::consecutive_alphabet_as_parallel_builder <std::uint8_t> builder;
			generate_alphabet(builder);
		}
	}
	
	
	void generate_context::calculate_segmentation_short_path(std::size_t const lb, std::size_t const rb)
	{
		segmentation_sp_context ctx(*this, lb, rb);
		ctx.process();
		std::cerr << " done.\nOutputting…";
		ctx.output_founders();
		
		if (m_segments_ostream_ptr)
			ctx.output_segments();
		
		std::cerr << std::endl;
		
		// Finish.
		cleanup();
		lb::log_time(std::cerr);
		std::cerr << "Done." << std::endl;
		exit(EXIT_SUCCESS);
	}
	
	
	void generate_context::calculate_segmentation_long_path(std::size_t const lb, std::size_t const rb)
	{
		assert(m_segment_length < rb - lb);
		assert(0 == lb); // FIXME: handle ranges that don't start from zero.
		
		auto *ctx(new segmentation_lp_context(*this, m_parallel_queue, m_serial_queue)); // Uses callbacks, deleted in the final one.
		ctx->generate_traceback(lb, rb);
	}
	
	
	void generate_context::check_traceback_size(segmentation_context &ctx)
	{
		if (! (ctx.max_segment_size() < sequence_count()))
		{
			// FIXME: cleanup
			std::cerr << "Unable to reduce the number of sequences; the maximum segment size is equal to the number of input sequences." << std::endl;
			exit(EXIT_FAILURE);
		}
	}
	
	
	void generate_context::context_did_finish_traceback(segmentation_sp_context &ctx)
	{
		check_traceback_size(ctx);
	}
	
	
	void generate_context::context_did_finish_traceback(segmentation_lp_context &ctx)
	{
		check_traceback_size(ctx);
		lb::log_time(std::cerr);
		std::cerr << "Updating the PBWT samples to traceback positions…" << std::endl;
		ctx.update_samples_to_traceback_positions();
	}
	
	
	void generate_context::context_did_update_pbwt_samples_to_traceback_positions(segmentation_lp_context &ctx)
	{
		segmentation_container container;
		
		lb::log_time(std::cerr);
		std::cerr << "Reducing the number of segments…" << std::endl;
		ctx.find_segments_greedy(container);
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
		assert(m_segmentation_istream.is_open());
		
		lb::log_time(std::cerr);
		std::cerr << "Loading the segmentation…" << std::endl;
		boost::archive::text_iarchive archive(m_segmentation_istream);
		archive >> m_alphabet;
		archive >> container;
	}

	
	void generate_context::save_segmentation_to_file(segmentation_container const &container)
	{
		assert(m_segmentation_ostream.is_open());
		
		lb::log_time(std::cerr);
		std::cerr << "Saving the segmentation…" << std::endl;
		boost::archive::text_oarchive archive(m_segmentation_ostream);
		archive << m_alphabet;
		archive << container;
	}
	
	
	void generate_context::join_segments_and_output(segmentation_container &&container)
	{
		lb::log_time(std::cerr);
		std::cerr << "Joining the remaining segments…" << std::endl;
		
		auto *ctx(new join_context(*this, m_parallel_queue, m_serial_queue, std::move(container), m_random_seed)); // Uses callbacks, deleted in the final one.
		ctx->join_segments_and_output(m_segment_joining_method);
	}
	
	
	void generate_context::context_did_output_founders(join_context &ctx)
	{
		if (m_segments_ostream_ptr)
			ctx.output_segments(m_segment_joining_method);
		
		ctx.cleanup();
		finish_lp();
	}
	
	
	void generate_context::finish_lp()
	{
		// Finish.
		cleanup();
		lb::log_time(std::cerr);
		std::cerr << "Done." << std::endl;
		exit(EXIT_SUCCESS);
	}
	
	
	void generate_context::calculate_segmentation(std::size_t const lb, std::size_t const rb)
	{
		lb::log_time(std::cerr);
		std::cerr << "Calculating the segmentation…" << std::flush;
		
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
		m_parallel_queue.reset(
			(
				m_use_single_thread ?
				dispatch_get_main_queue() :
				dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0)
			),
			true
		);
		
		m_serial_queue.reset(
			(
				m_use_single_thread ?
				dispatch_get_main_queue() :
				dispatch_queue_create("fi.iki.tsnorri.processing-queue", DISPATCH_QUEUE_SERIAL)
			),
			m_use_single_thread
		);
		
		if (segmentation_input_path)
			lb::open_file_for_reading(segmentation_input_path, m_segmentation_istream);
		else if (segmentation_output_path)
			lb::open_file_for_writing(segmentation_output_path, m_segmentation_ostream, lb::writing_open_mode::CREATE);
		
		if (output_founders_path)
			lb::open_file_for_writing(output_founders_path, m_founders_ostream, lb::writing_open_mode::CREATE);
		
		if (output_segments_path)
		{
			if ('-' == output_segments_path[0] && '\0' == output_segments_path[1])
				m_segments_ostream_ptr = &std::cerr;
			else
			{
				lb::open_file_for_writing(output_segments_path, m_segments_ostream, lb::writing_open_mode::CREATE);
				m_segments_ostream_ptr = &m_segments_ostream;
			}
		}
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
			generate_alphabet();
			
			auto const sequence_length(m_sequences.front().size());
			if (0 == m_pbwt_sample_rate)
				m_pbwt_sample_rate = std::ceil(std::sqrt(sequence_length));
			
			calculate_segmentation(0, sequence_length);
		}
	}
}
