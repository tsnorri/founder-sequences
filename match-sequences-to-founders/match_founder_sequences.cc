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
#include <libbio/vector_source.hh>
#include <mutex>
#include <numeric>
#include <string>
#include <vector>


namespace fseq	= founder_sequences;
namespace lb	= libbio;


namespace {
	
	template <typename t_vector_source>
	class line_reader_cb
	{
	protected:
		typedef t_vector_source							vector_source;
		typedef typename vector_source::vector_type		vector_type;
		
	protected:
		std::unique_ptr <vector_type>	m_sequence_ptr;
		
	public:
		std::unique_ptr <vector_type> &sequence_ptr() { return m_sequence_ptr; }
		
		void handle_sequence(
			uint32_t line,
			std::unique_ptr <vector_type> &sequence_ptr,
			std::size_t const &seq_length,
			vector_source &vector_source
		)
		{
			m_sequence_ptr = std::move(sequence_ptr);
			m_sequence_ptr->resize(seq_length);
			sequence_ptr.reset();
			vector_source.set_vector_length(seq_length);
			// Don't return the vector to the source.
		}
		
		void start() {}
		void finish() {}
	};
	
	
	class match_context final
	{
	protected:
		typedef lb::vector_source <std::vector <std::uint8_t>>	vector_source;
		
	protected:
		lb::dispatch_ptr <dispatch_queue_t>		m_queue;
		lb::dispatch_ptr <dispatch_group_t>		m_matching_group;
		lb::dispatch_ptr <dispatch_semaphore_t>	m_reading_semaphore;
		std::mutex								m_output_mutex;
		
		vector_source							m_vector_source;
		
		fseq::sequence_vector					m_founder_matrix;
		std::vector <std::string>				m_sequence_paths;
		
		bool									m_use_single_thread{false};
	
	public:
		match_context() = delete;
		match_context(
			std::vector <std::string> &&sequence_paths,
			fseq::sequence_vector &&founder_matrix,
			bool use_single_thread
		):
			m_founder_matrix(std::move(founder_matrix)),
			m_sequence_paths(std::move(sequence_paths)),
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
			auto const &founder_seq(m_founder_matrix[seq_idx]);
			auto const founder_char(founder_seq[char_idx]);
			if (founder_char == c)
			{
				++retval;
				dst_seq_indices.push_back(seq_idx);
			}
		}
		
		return retval;
	}
	
	
	void match_context::match_sequence_and_report(std::vector <uint8_t> const &sequence, std::size_t const seq_idx)
	{
		// Create an initial permutation of all the founder sequences. Then filter those sequences
		// while iterating the characters in the given original sequence.
		std::vector <std::size_t> seq_indices(m_founder_matrix.size());
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
				{
					std::lock_guard <std::mutex> lock(m_output_mutex);
					std::cout << seq_idx << lb << '\t' << chr_idx << '\t';
					std::copy(
						seq_indices.cbegin(),
						seq_indices.cend(),
						std::experimental::make_ostream_joiner(std::cout, ",")
					);
					std::cout << '\n';
				}
				
				// Re-initialize for the new range.
				lb = chr_idx;
				count = m_founder_matrix.size();
				seq_indices.resize(count);
				std::iota(seq_indices.begin(), seq_indices.end(), 0);
				
				// Re-match with the new range.
				dst_count = compare_to_sequences(seq_indices, c, chr_idx, dst_seq_indices);
				if (0 == dst_count)
				{
					std::lock_guard <std::mutex> lock(m_output_mutex);
					std::cerr << "Error: character '" << c << "' (" << +c << ") at " << seq_idx << ':' << chr_idx << " not found in the founders." << std::endl;
				}
			}
			
			using std::swap;
			swap(count, dst_count);
			swap(seq_indices, dst_seq_indices);
			++chr_idx;
		}
		
		if (0 == count)
		{
			auto const c(sequence[chr_idx - 1]);
			
			std::lock_guard <std::mutex> lock(m_output_mutex);
			std::cerr << "Error: character '" << c << "' (" << +c << ") at " << seq_idx << ':' << chr_idx << " not found in the founders." << std::endl;
		}
		else
		{
			// Output the current range.
			std::lock_guard <std::mutex> lock(m_output_mutex);
			std::cout << seq_idx << lb << '\t' << chr_idx << '\t';
			std::copy(
				seq_indices.cbegin(),
				seq_indices.cend(),
				std::experimental::make_ostream_joiner(std::cout, ",")
			);
			std::cout << '\n';
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
					
					typedef line_reader_cb <vector_source> line_reader_cb;
					typedef lb::line_reader <vector_source, line_reader_cb, 0> line_reader;
				
					lb::file_istream input_stream;
					lb::open_file_for_reading(path.c_str(), input_stream);
				
					line_reader reader;
					line_reader_cb cb;
					reader.read_from_stream(input_stream, m_vector_source, cb);
					
					using std::swap;
					swap(sequence_ptr, cb.sequence_ptr());
					
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
	
	
	void read_file_contents(char const *path, std::vector <std::string> lines)
	{
		lb::file_istream stream;
		lb::open_file_for_reading(path, stream);
		
		std::string line;
		while (std::getline(stream, line))
			lines.push_back(line);
	}
	
	
	bool read_single_line(char const *path, std::string &line)
	{
		lb::file_istream stream;
		lb::open_file_for_reading(path, stream);
		return bool(std::getline(stream, line));
	}
	
	
	void read_input(
		char const *sequences_list_path,
		char const *founders_path,
		std::vector <std::string> &sequence_file_paths,
		fseq::sequence_vector &founder_matrix
	)
	{
		std::cerr << "Reading sequence paths…" << std::endl;
		read_file_contents(sequences_list_path, sequence_file_paths);
		
		std::cerr << "Reading founders…" << std::endl;
		{
			typedef lb::vector_source <std::vector <std::uint8_t>> vector_source;
			typedef fseq::line_reader_cb <vector_source> line_reader_cb;
			typedef lb::line_reader <vector_source, line_reader_cb, 0> line_reader;
		
			vector_source vs;
			line_reader reader;
			line_reader_cb cb(founder_matrix);
		
			lb::file_istream stream;
			lb::open_file_for_reading(founders_path, stream);
			reader.read_from_stream(stream, vs, cb);
		}
	}
}


namespace founder_sequences {
	void match_founder_sequences(
		char const *sequences_list_path,
		char const *founders_path,
		bool const single_threaded
	)
	{
		std::vector <std::string> sequence_file_paths;
		fseq::sequence_vector founder_matrix;
		read_input(sequences_list_path, founders_path, sequence_file_paths, founder_matrix);
		
		std::cerr << "Matching founders with sequences…" << std::endl;
		auto *ctx(new match_context(std::move(sequence_file_paths), std::move(founder_matrix), single_threaded));
	
		ctx->prepare();
		ctx->match();
	}
}
