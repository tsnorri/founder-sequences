/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <algorithm>
#include <boost/range/combine.hpp>
#include <experimental/iterator>
#include <founder_sequences/founder_sequences.hh>
#include <libbio/dispatch_fn.hh>
#include <libbio/file_handling.hh>
#include <libbio/line_reader.hh>
#include <libbio/sequence_reader/sequence_reader.hh>
#include <libbio/vector_source.hh>
#include <mutex>
#include <numeric>
#include <string>
#include <vector>


namespace fseq	= founder_sequences;
namespace lb	= libbio;
namespace lsr	= libbio::sequence_reader;


namespace {
	
	class match_context final
	{
	protected:
		typedef lb::vector_source <std::vector <std::uint8_t>>	vector_source;
		
	protected:
		lb::dispatch_ptr <dispatch_queue_t>			m_queue;
		lb::dispatch_ptr <dispatch_group_t>			m_matching_group;
		lb::dispatch_ptr <dispatch_semaphore_t>		m_reading_semaphore;
		std::mutex									m_output_mutex;
		
		vector_source								m_vector_source;
		
		std::vector <std::string>					m_sequence_paths;
		std::unique_ptr <lsr::sequence_container>	m_founder_container;
		lsr::sequence_vector						m_founders;
		
		bool										m_use_single_thread{false};
	
	public:
		match_context() = delete;
		match_context(
			std::vector <std::string> &&sequence_paths,
			std::unique_ptr <lsr::sequence_container> &&founder_container,
			bool use_single_thread
		):
			m_sequence_paths(std::move(sequence_paths)),
			m_founder_container(std::move(founder_container)),
			m_use_single_thread(use_single_thread)
		{
		}
		
		void prepare();
		void match();
		void cleanup() { delete this; }
		
	protected:
		std::size_t compare_to_sequences(
			std::vector <std::size_t> const &seq_indices,
			char const c,
			std::size_t const char_idx,
			std::vector <std::size_t> &dst_seq_indices
		) const;
		
		void match_sequence_and_report(std::vector <uint8_t> const &sequence, std::size_t const seq_idx);
		void output_range(std::size_t const seq_idx, std::size_t const lb, std::size_t const chr_idx, std::vector <std::size_t> const &seq_indices);
		void output_character_not_found_error(char const c, std::size_t const seq_idx, std::size_t const chr_idx);
	};
	
	
	void match_context::prepare()
	{
		m_queue.reset(
			(
				m_use_single_thread ?
				dispatch_get_main_queue() :
				dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0)
			),
			true
		);
		
		m_matching_group.reset(dispatch_group_create());
		m_reading_semaphore.reset(dispatch_semaphore_create(1));
		
		m_founder_container->to_spans(m_founders);
	}
	
	
	std::size_t match_context::compare_to_sequences(
		std::vector <std::size_t> const &seq_indices,
		char const c,
		std::size_t const char_idx,
		std::vector <std::size_t> &dst_seq_indices
	) const
	{
		std::size_t retval(0);
		dst_seq_indices.clear(); // Leaves the capacity() of std::vector unchanged.
		
		for (auto const seq_idx : seq_indices)
		{
			auto const &founder_seq(m_founders[seq_idx]);
			auto const founder_char(founder_seq[char_idx]);
			if (founder_char == c)
			{
				++retval;
				dst_seq_indices.push_back(seq_idx);
			}
		}
		
		return retval;
	}
	
	
	void match_context::output_range(std::size_t const seq_idx, std::size_t const lb, std::size_t const chr_idx, std::vector <std::size_t> const &seq_indices)
	{
		std::lock_guard <std::mutex> lock(m_output_mutex);
		std::cout << seq_idx << '\t' << lb << '\t' << chr_idx << '\t';
		std::copy(
			seq_indices.cbegin(),
			seq_indices.cend(),
			std::experimental::make_ostream_joiner(std::cout, ",")
		);
		std::cout << '\n';
	}
	
	
	void match_context::output_character_not_found_error(char const c, std::size_t const seq_idx, std::size_t const chr_idx)
	{
		std::lock_guard <std::mutex> lock(m_output_mutex);
		std::cerr << "Error: character '" << c << "' (" << +c << ") at " << seq_idx << ':' << chr_idx << " not found in the founders." << std::endl;
	}

	
	void match_context::match_sequence_and_report(std::vector <uint8_t> const &sequence, std::size_t const seq_idx)
	{
		// Create an initial permutation of all the founder sequences. Then filter those sequences
		// while iterating the characters in the given original sequence.
		std::vector <std::size_t> seq_indices(m_founders.size());
		std::vector <std::size_t> dst_seq_indices;
		std::iota(seq_indices.begin(), seq_indices.end(), 0);
		
		std::size_t lb(0);
		std::size_t count(seq_indices.size());
		std::size_t dst_count(0);
		std::size_t chr_idx(0);
		for (auto const c : sequence)
		{
			dst_count = compare_to_sequences(seq_indices, c, chr_idx, dst_seq_indices);
			if (0 == dst_count)
			{
				// Output the current range.
				output_range(seq_idx, lb, chr_idx, seq_indices);
				
				// Re-initialize for the new range.
				lb = chr_idx;
				count = m_founders.size();
				seq_indices.resize(count);
				std::iota(seq_indices.begin(), seq_indices.end(), 0);
				
				// Re-match with the new range.
				dst_count = compare_to_sequences(seq_indices, c, chr_idx, dst_seq_indices);
				if (0 == dst_count)
					output_character_not_found_error(c, seq_idx, chr_idx);
			}
			
			using std::swap;
			swap(count, dst_count);
			swap(seq_indices, dst_seq_indices);
			++chr_idx;
		}
		
		if (0 == count)
		{
			auto const c(sequence[chr_idx - 1]);
			output_character_not_found_error(c, seq_idx, chr_idx);
		}
		else
		{
			// Output the current range.
			output_range(seq_idx, lb, chr_idx, seq_indices);
		}
	}
	
	
	void match_context::match()
	{
		std::cout << "SEQUENCE_INDEX" "\t" "LB" "\t" "RB" "\t" "FOUNDER_INDICES" "\n";
		
		std::size_t seq_idx(0);
		for (auto const &path : m_sequence_paths)
		{
			auto const current_seq_idx(seq_idx++);
			dispatch_group_async(*m_matching_group, *m_queue, ^{
				
				std::unique_ptr <typename vector_source::vector_type> sequence_ptr;
				
				{
					// Read one original sequence from the disk.
					// Make just one thread read from the disk at a time.
					dispatch_semaphore_wait(*m_reading_semaphore, DISPATCH_TIME_FOREVER);
					
					m_vector_source.get_vector(sequence_ptr);
					lb::file_istream input_stream;
					lb::open_file_for_reading(path.c_str(), input_stream);
					lb::read_from_stream(input_stream, *sequence_ptr);
					
					dispatch_semaphore_signal(*m_reading_semaphore);
				}
				
				// Handle the sequence.
				match_sequence_and_report(*sequence_ptr, current_seq_idx);
				
				// Make the buffer available for the next thread.
				m_vector_source.put_vector(sequence_ptr);
			});
		}
		
		// Output each founder sequence by following the matched pairs.
		dispatch_group_notify(*m_matching_group, dispatch_get_main_queue(), ^{
			std::cout << std::flush;
			
			// Exit.
			cleanup();
			exit(EXIT_SUCCESS);
		});
	}
}


namespace founder_sequences {
	void match_founder_sequences(
		char const *sequences_path,
		char const *founders_path,
		lsr::input_format const founders_format,
		bool const single_threaded
	)
	{
		std::vector <std::string> paths;
		std::unique_ptr <lsr::sequence_container> founder_container;
		
		std::cerr << "Reading sequence paths…" << std::endl;
		lsr::read_list_file(sequences_path, paths);
		std::cerr << "Reading founders…" << std::endl;
		lsr::read_input(founders_path, founders_format, founder_container);
		
		std::cerr << "Matching founders with sequences…" << std::endl;
		auto *ctx(new match_context(std::move(paths), std::move(founder_container), single_threaded));
		ctx->prepare();
		ctx->match();
	}
}
