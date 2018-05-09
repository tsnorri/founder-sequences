/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <algorithm>
#include <boost/range/combine.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/range/combine.hpp>
#include <experimental/iterator>
#include <founder_sequences/founder_sequences.hh>
#include <founder_sequences/merge_segments_task.hh>
#include <founder_sequences/rmq.hh>
#include <libbio/algorithm.hh>
#include <libbio/consecutive_alphabet.hh>
#include <libbio/counting_sort.hh>
#include <libbio/dispatch_fn.hh>
#include <libbio/fasta_reader.hh>
#include <libbio/file_handling.hh>
#include <libbio/line_reader.hh>
#include <libbio/pbwt.hh>
#include <libbio/vector_source.hh>
#include <memory>
#include <sdsl/int_vector.hpp>
#include <sdsl/rmq_support.hpp>
#include <vector>

#define PRINT_SEGMENTATION_TRACEBACK	0
#define PRINT_SEGMENT_TEXTS				0


namespace fseq	= founder_sequences;
namespace lb	= libbio;


namespace {
	
	struct index_divergence_pair;
	
	typedef lb::consecutive_alphabet_as <std::uint8_t>	alphabet_type;
	typedef std::vector <index_divergence_pair>			index_divergence_pair_vector;
	
	typedef lb::pbwt::pbwt_context <
		sdsl::int_vector <32>,
		sdsl::int_vector <32>,
		sdsl::range_maximum_sct <>::type,
		sdsl::int_vector <32>,
		fseq::sequence_vector,
		alphabet_type,
		uint32_t
	> pbwt_context;
	
	
	struct segmentation_dp_arg
	{
		std::size_t		lb{0};	// Inclusive.
		std::size_t 	rb{0};	// Exclusive, i.e. [ ) style range.
		std::uint32_t	segment_max_size{UINT32_MAX}; // Maximum of the current segment and the previous one.
		std::uint32_t	segment_size{UINT32_MAX};

		segmentation_dp_arg() = default;
		
		segmentation_dp_arg(
			std::size_t const lb_,
			std::size_t const rb_,
			std::uint32_t const segment_max_size_,
			std::uint32_t const segment_size_
		):
			lb(lb_),
			rb(rb_),
			segment_max_size(segment_max_size_),
			segment_size(segment_size_)
		{
			assert(lb <= rb);
		}
		
		bool operator<(segmentation_dp_arg const &rhs) const { return segment_max_size < rhs.segment_max_size; }
		std::size_t text_length() const { return rb - lb; }
	};
	
	std::ostream &operator<<(std::ostream &os, segmentation_dp_arg const &dp_arg);
	
	
	class generate_context final
	{
	protected:
		typedef std::vector <segmentation_dp_arg>							segmentation_traceback_vector;
		typedef fseq::rmq <segmentation_traceback_vector, std::less <>, 64>	segmentation_traceback_vector_rmq;
		
	protected:
		lb::dispatch_ptr <dispatch_queue_t>							m_queue;
		lb::dispatch_ptr <dispatch_group_t>							m_matching_group;
		
		fseq::sequence_vector										m_sequences;
		alphabet_type												m_alphabet;
		segmentation_traceback_vector								m_segmentation_traceback;
		fseq::segment_text_matrix									m_segment_texts;
		
		std::vector <std::unique_ptr <fseq::merge_segments_task>>	m_merge_tasks;
		std::vector <fseq::matching_vector>							m_matchings;
		
		std::size_t													m_segment_length{0};
		fseq::segment_joining										m_segment_joining_method{fseq::segment_joining::MATCHING};
		bool														m_use_single_thread{false};
		
	public:
		generate_context(
			std::size_t const segment_length,
			fseq::segment_joining const segment_joining_method,
			bool const use_single_thread
		):
			m_segment_length(segment_length),
			m_segment_joining_method(segment_joining_method),
			m_use_single_thread(use_single_thread)
		{
		}
		
		generate_context(generate_context const &) = delete;
		generate_context(generate_context &&) = delete;
		
		void prepare();
		void load_and_generate(
			char const *input_path,
			fseq::input_format const input_file_format,
			char const *output_segments_path
		);
		void cleanup() { delete this; }
		
		// For debugging.
		void print_segment_texts() const;
		void print_segmentation_traceback(std::size_t const idx, segmentation_dp_arg const &arg) const;
		
	protected:
		void load_input(char const *input_path, fseq::input_format const input_file_format);
		void load_input_fasta(lb::file_istream &stream);
		void load_input_list_file(lb::file_istream &stream);
		void check_input() const;
		void generate_alphabet();
		void generate_founders(std::size_t const lb, std::size_t const rb);
		
		void calculate_segmentation(std::size_t const lb, std::size_t const rb);
		void calculate_segmentation_short_path(std::size_t const lb, std::size_t const rb);
		void calculate_segmentation_long_path(std::size_t const lb, std::size_t const rb);
		
		void calculate_segmentation_lp_dp_arg(
			pbwt_context const &pbwt_ctx,
			segmentation_traceback_vector const &segmentation_traceback_dp,
			segmentation_traceback_vector_rmq const &segmentation_traceback_rmq,
			std::size_t const lb,	// Inclusive.
			std::size_t const text_pos,
			segmentation_dp_arg &min_arg
		);
		
		void calculate_segmentation_lp_fill_traceback(
			std::size_t const lb,
			std::size_t const rb,
			pbwt_context &pbwt_ctx,
			segmentation_traceback_vector &segmentation_traceback_res
		);
		
		void calculate_segmentation_lp_create_segments(
			segmentation_traceback_vector const &segmentation_traceback,
			std::size_t const max_segment_size,
			pbwt_context &pbwt_ctx,
			fseq::segment_text_matrix &segment_texts
		);
		
		void calculate_segmentation_lp_copy_segments(
			std::size_t const max_segment_size,
			fseq::segment_text_matrix &segment_texts,
			std::size_t const sequence_count
		);
		
		void fill_dp_arg(
			std::size_t const text_pos,
			std::size_t const dp_lb,
			std::size_t const dp_rb,
			uint32_t const rhs,
			segmentation_traceback_vector const &segmentation_traceback,
			segmentation_traceback_vector_rmq const &segmentation_traceback_rmq,
			segmentation_dp_arg &dst_arg
		) const;
		
		void copy_segment_texts(
			pbwt_context const &pbwt_ctx,
			std::size_t const lb,
			std::size_t const rb,								// [ ) style range.
			fseq::segment_text_vector &current_segment_texts,	// Buffer provided by caller.
			std::vector <std::size_t> &suffix_numbers,			// Buffer provided by caller.
			fseq::segment_text_matrix &segment_texts
		);
		
		void output_segments(lb::file_ostream &stream) const;
		void output_founders_short_path() const;
		void match_segments_random_and_output();
		void match_segments_sd_and_output();
	};
	
	
	void generate_context::load_input_list_file(lb::file_istream &list_stream)
	{
		// Prepare for reading the sequences.
		typedef lb::vector_source <std::vector <std::uint8_t>> vector_source;
		typedef fseq::line_reader_cb <vector_source> line_reader_cb;
		typedef lb::line_reader <vector_source, line_reader_cb, 0> line_reader;
		
		vector_source vs;
		line_reader reader;
		line_reader_cb cb(m_sequences);
		
		// Read the input file names and handle each file.
		std::string path;
		while (std::getline(list_stream, path))
		{
			lb::file_istream stream;
			lb::open_file_for_reading(path.c_str(), stream);
			reader.read_from_stream(stream, vs, cb);
		}
	}
	
	
	void generate_context::load_input_fasta(lb::file_istream &stream)
	{
		// Read the sequences.
		typedef lb::vector_source <std::vector <std::uint8_t>> vector_source;
		typedef fseq::fasta_reader_cb <vector_source> fasta_reader_cb;
		typedef lb::fasta_reader <vector_source, fasta_reader_cb, 0> fasta_reader;
		
		vector_source vs;
		fasta_reader reader;
		fasta_reader_cb cb(m_sequences);
		reader.read_from_stream(stream, vs, cb);
	}
	
	
	void generate_context::load_input(char const *input_path, fseq::input_format const input_file_format)
	{
		std::cerr << "Loading the input…" << std::endl;
		
		// Open the input file.
		lb::file_istream input_stream;
		lb::open_file_for_reading(input_path, input_stream);
		
		switch (input_file_format)
		{
			case fseq::input_format::FASTA:
				load_input_fasta(input_stream);
				break;
			
			case fseq::input_format::LIST_FILE:
				load_input_list_file(input_stream);
				break;
			
			default:
				lb::fail("Unexpected input file format.");
				break;
		}
		
		if (0 == m_sequences.size())
		{
			std::cerr << "The input file contained no sequences." << std::endl;
			exit(EXIT_SUCCESS);
		}
	}
	
	
	std::ostream &operator<<(std::ostream &os, segmentation_dp_arg const &dp_arg)
	{
		os << '[' << dp_arg.lb << ", " << dp_arg.rb << ") segment_max_size: " << dp_arg.segment_max_size << " segment_size: " << dp_arg.segment_size;
		return os;
	}
	
	
	void generate_context::check_input() const
	{
		std::cerr << "Checking the input…" << std::endl;

		// Check that all the vectors have equal lengths.
		auto const sequence_length(m_sequences.front().size());
		bool stop(false);
		std::size_t i(0);
		for (auto const &vec : m_sequences)
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
		std::cerr << "Generating a compressed alphabet…" << std::endl;
		for (auto const &vec : m_sequences)
			m_alphabet.prepare(vec);
		m_alphabet.compress();
	}
	
	
	void generate_context::fill_dp_arg(
		std::size_t const text_pos,
		std::size_t const dp_lb,
		std::size_t const dp_rb,
		std::uint32_t const rhs,
		segmentation_traceback_vector const &segmentation_traceback_dp,
		segmentation_traceback_vector_rmq const &segmentation_traceback_rmq,
		segmentation_dp_arg &dst_arg
	) const
	{
		// Convert to segmentation_traceback indexing.
		auto const dp_lb_c(1 + dp_lb - m_segment_length);
		auto const dp_rb_c(1 + dp_rb - m_segment_length);
		auto const idx(segmentation_traceback_rmq(dp_lb_c, dp_rb_c));
		
		auto const &boundary_segment(segmentation_traceback_dp[idx]);
		auto const lhs(boundary_segment.segment_max_size);	// Corresponds to M(h).
		
		segmentation_dp_arg arg(idx + m_segment_length, 1 + text_pos, std::max(lhs, rhs), rhs);
		dst_arg = std::move(arg);
	}
	
	
	void generate_context::copy_segment_texts(
		pbwt_context const &pbwt_ctx,
		std::size_t const lb,
		std::size_t const rb,								// [ ) style range.
		fseq::segment_text_vector &current_segment_texts,	// Buffer provided by caller.
		std::vector <std::size_t> &suffix_numbers,			// Buffer provided by caller.
		fseq::segment_text_matrix &segment_texts
	
	)
	{
		// Get the text segments.
		{
			current_segment_texts.clear();
			std::size_t suffix_idx(0);
			//pbwt_ctx.print_vectors();
			for (auto const &zipped : boost::combine(pbwt_ctx.input_permutation(), pbwt_ctx.input_divergence()))
			{
				auto const sequence_idx(zipped.get <0>());
				auto const divergence(zipped.get <1>());
				
				if (! (divergence <= lb))
				{
					// Output the suffix since it is distinct from the previous one.
					// Number the suffixes as they are being output.
					auto const &sequence(m_sequences[sequence_idx]);
					fseq::segment_text seg_text;
					
					{
						std::string substring(sequence.cbegin() + lb, sequence.cbegin() + rb);
						seg_text.text = std::move(substring);
					}
					
					current_segment_texts.emplace_back(std::move(seg_text));
					++suffix_idx;
				}
				
				// Memorize the suffix number of the current string.
				assert(sequence_idx < suffix_numbers.size());
				suffix_numbers[sequence_idx] = suffix_idx - 1;
			}
		}
		
		// Annotate the segment with the sequence numbers in sorted order.
		{
			auto const sequence_count(pbwt_ctx.size());
			for (std::size_t i(0); i < sequence_count; ++i)
			{
				auto const suffix_idx(suffix_numbers[i]);
				assert(suffix_idx < current_segment_texts.size());
				current_segment_texts[suffix_idx].sequence_indices.push_back(i);
			}
		}
		
		// Store the current segments.
		segment_texts.emplace_back(std::move(current_segment_texts));
	}
	
	
	void generate_context::calculate_segmentation_short_path(std::size_t const lb, std::size_t const rb)
	{
		pbwt_context pbwt_ctx(m_sequences, m_alphabet);
		pbwt_ctx.prepare();
		
		for (std::size_t j(lb); j < rb; ++j)
		{
			pbwt_ctx.build_prefix_and_divergence_arrays(j);
			pbwt_ctx.update_divergence_value_counts();
			pbwt_ctx.swap_input_and_output();
		}
		
		auto const seq_count(pbwt_ctx.size());
		std::size_t segment_size_diff(0);
		auto const &counts(pbwt_ctx.divergence_value_counts());
		auto const begin(counts.cbegin_pairs());
		if (counts.cend_pairs() != begin && lb == begin->first)
			segment_size_diff = begin->second;
		
		auto const partition_size(seq_count - segment_size_diff);
		m_segmentation_traceback.emplace_back(lb, rb, partition_size, partition_size);
		
		// Read the texts.
		fseq::segment_text_vector current_segment_texts;
		std::vector <std::size_t> suffix_numbers(pbwt_ctx.size(), 0); // Suffix numbers by string index.
		copy_segment_texts(
			pbwt_ctx,
			lb,
			rb,
			current_segment_texts,
			suffix_numbers,
			m_segment_texts
		);
	}
	
	
	void generate_context::calculate_segmentation_lp_dp_arg(
		pbwt_context const &pbwt_ctx,
		segmentation_traceback_vector const &segmentation_traceback_dp,
		segmentation_traceback_vector_rmq const &segmentation_traceback_rmq,
		std::size_t const lb,	// Inclusive.
		std::size_t const text_pos,
		segmentation_dp_arg &min_arg
	)
	{
		auto const seq_count(pbwt_ctx.size());
		
		// Take the current divergence values and their counts in divergence value order.
		auto it(pbwt_ctx.divergence_value_counts().cbegin_pairs());
		auto const end(pbwt_ctx.divergence_value_counts().cend_pairs());
		assert(it != end);

		// Find the minimum, similar to Ukkonen's equation 1.
		std::size_t segment_size_diff(it->second);

		std::size_t dp_lb(0);
		std::size_t dp_rb(it->first); // Initial value not used.
		
		//pbwt_ctx.print_vectors();
		
		// Consider the whole range in case the segment size is smaller than seq_count.
		if (lb == dp_rb)
		{
			auto const segment_size(seq_count - segment_size_diff);
			segmentation_dp_arg const current_arg(lb, 1 + text_pos, segment_size, segment_size);
			if (current_arg < min_arg)
				min_arg = current_arg;
			
			++it;
			dp_rb = it->first;
			segment_size_diff += it->second;

		}
		
		++it;
		while (true)
		{
			if (it == end)
				break;
			
			// Range of possible cutting points.
			dp_lb = dp_rb - 1;
			dp_rb = it->first - 1;
			std::size_t dp_rb_c(dp_rb);
			assert(0 != it->first);
			assert(dp_lb < dp_rb);
			assert(dp_rb <= 1 + text_pos);
			
			// Check if the range needs to and can be contracted.
			// First, verify that the current segment fits after the DP search range.
			if (1 + text_pos - m_segment_length < dp_rb_c)
				dp_rb_c = 1 + text_pos - m_segment_length;
			
			// Second, verify that one segment fits before the DP search range.
			if (dp_lb < lb + m_segment_length - 1)
			{
				if (lb + m_segment_length - 1 < dp_rb_c)
					dp_lb = lb + m_segment_length - 1;
				else
					goto continue_loop;
			}
			
			// Check that the range is still valid.
			if (dp_lb < dp_rb_c)
			{
				// Convert to segmentation_traceback indexing.
				auto const dp_lb_tb(1 + dp_lb - m_segment_length);
				auto const dp_rb_tb(1 + dp_rb_c - m_segment_length);
				auto const idx(segmentation_traceback_rmq(dp_lb_tb, dp_rb_tb));
				
				auto const &boundary_segment(segmentation_traceback_dp[idx]);
				auto const lhs(boundary_segment.segment_max_size);
				auto const rhs(seq_count - segment_size_diff);
				
				segmentation_dp_arg current_arg(idx + m_segment_length, 1 + text_pos, lb::max_ct(lhs, rhs), rhs);
				if (current_arg < min_arg)
					min_arg = current_arg;
			}
			
			// Next iteration.
		continue_loop:
			segment_size_diff += it->second;
			++it;
			
			//pbwt_ctx.print_vectors();
		}
	}
	
	
	void generate_context::calculate_segmentation_lp_fill_traceback(
		std::size_t const lb,	// Inclusive.
		std::size_t const rb,	// Exclusive, i.e. [ ) .
		pbwt_context &pbwt_ctx,
		segmentation_traceback_vector &segmentation_traceback_res
	)
	{
		pbwt_ctx.prepare();
		
		auto const seq_length(m_sequences.front().size());
		auto const seq_count(pbwt_ctx.size());
		auto const dp_size(seq_length - m_segment_length + 1);
		
		segmentation_traceback_vector segmentation_traceback_dp(dp_size);	// Values shifted to the left by m_segment_length (L) since the first L columns have the same value anyway.
		segmentation_traceback_vector_rmq segmentation_traceback_rmq(segmentation_traceback_dp);
		
		// Calculate the first L - 1 columns, which gives the required result for calculating M(L).
		// j is 0-based, m_segment_length is 1-based.
		std::size_t j(lb);
		
		while (j < m_segment_length - 1)
		{
			pbwt_ctx.build_prefix_and_divergence_arrays(j);
			pbwt_ctx.update_divergence_value_counts();
			pbwt_ctx.swap_input_and_output();
			++j;
		}
		
		//pbwt_ctx.print_vectors();
		
		// Calculate the columns L to 2L - 1 (1-based), which gives the size of one segment for each column since the minimum
		// length of two segments is 2L.
		for (std::size_t const limit(std::min(2 * m_segment_length, rb - m_segment_length) - 1); j < limit; ++j)
		{
			pbwt_ctx.build_prefix_and_divergence_arrays(j);
			pbwt_ctx.update_divergence_value_counts();
			
			// Calculate the segment size by finding the range of the relevant key
			// (which is “0” in this case b.c. the column count is less than 2L,
			// so the range is at most [counts.begin(), counts.begin() + 1))
			// and subtracting the count from the sequence count.
			std::size_t segment_size_diff(0);
			auto const &counts(pbwt_ctx.divergence_value_counts());
			auto const begin(counts.cbegin_pairs());
			if (counts.cend_pairs() != begin && 0 == begin->first)
				segment_size_diff = begin->second;
			
			auto const idx(j + 1 - m_segment_length);
			auto const segment_size(seq_count - segment_size_diff);
			segmentation_dp_arg const current_arg(lb, 1 + j, segment_size, segment_size);
			print_segmentation_traceback(idx, current_arg);
			segmentation_traceback_dp[idx] = current_arg;
			segmentation_traceback_rmq.update(idx);
			
			pbwt_ctx.swap_input_and_output();
		}
		
		// Calculate the columns (L or) 2L to n (1-based). This results in multiple segments the minimum length of each is L.
		// Both rb an m_segment_length are 1-based.
		for (std::size_t const limit(rb - m_segment_length); j < limit; ++j)
		{
			pbwt_ctx.build_prefix_and_divergence_arrays(j);
			pbwt_ctx.update_divergence_value_counts();
			
			// Use the texts up to this point as the initial value.
			segmentation_dp_arg min_arg(lb, 1 + j, seq_count, seq_count);
			calculate_segmentation_lp_dp_arg(pbwt_ctx, segmentation_traceback_dp, segmentation_traceback_rmq, lb, j, min_arg);
			
			auto const idx(j + 1 - m_segment_length);
			print_segmentation_traceback(idx, min_arg);
			segmentation_traceback_dp[idx] = min_arg;
			segmentation_traceback_rmq.update(idx);
			
			pbwt_ctx.swap_input_and_output();
		}
		
		// Calculate the size of the final segment.
		{
			while (j < rb)
			{
				pbwt_ctx.build_prefix_and_divergence_arrays(j);
				pbwt_ctx.update_divergence_value_counts();
				pbwt_ctx.swap_input_and_output();
				++j;
			}
			
			//pbwt_ctx.print_vectors();
			
			segmentation_dp_arg min_arg(lb, j, seq_count, seq_count);
			calculate_segmentation_lp_dp_arg(pbwt_ctx, segmentation_traceback_dp, segmentation_traceback_rmq, lb, j - 1, min_arg);
			
			// Position of the last traceback argument.
			auto const idx(segmentation_traceback_dp.size() - 1);
			assert(rb - m_segment_length == idx);

			print_segmentation_traceback(idx, min_arg);
			segmentation_traceback_dp[idx] = min_arg;
		}
		
		// Follow the traceback.
		{
			assert(segmentation_traceback_dp.size());
			segmentation_traceback_res.clear();
			std::size_t arg_idx(segmentation_traceback_dp.size() - 1);
			while (true)
			{
				auto const &current_arg(segmentation_traceback_dp[arg_idx]);
				// Fill in reverse order.
				segmentation_traceback_res.push_back(current_arg);
				
				auto const next_pos(current_arg.lb);
				if (lb == next_pos)
					break;
				
				assert(m_segment_length <= next_pos);
				arg_idx = next_pos - m_segment_length;
			}
			
			// Reverse the filled traceback.
			std::reverse(segmentation_traceback_res.begin(), segmentation_traceback_res.end());
		}
	}
	
	
	void generate_context::calculate_segmentation_lp_create_segments(
		segmentation_traceback_vector const &segmentation_traceback,
		std::size_t const max_segment_size,
		pbwt_context &pbwt_ctx,
		fseq::segment_text_matrix &segment_texts
	)
	{
		// Create the segments, duplicate the items in proportion when the segment size is smaller than the maximum.
		fseq::segment_text_vector current_segment_texts;
		std::vector <std::size_t> suffix_numbers(pbwt_ctx.size(), 0); // Suffix numbers by string index.
		std::size_t j(0);
		
		pbwt_ctx.prepare();
		//pbwt_ctx.print_vectors();
		for (auto const &traceback_arg : segmentation_traceback)
		{
			// Re-calculate the PBWT up to the current column.
			while (j < traceback_arg.rb)
			{
				pbwt_ctx.build_prefix_and_divergence_arrays(j);
				pbwt_ctx.swap_input_and_output();
				++j;
			}
			
			//pbwt_ctx.print_vectors();
			copy_segment_texts(
				pbwt_ctx,
				traceback_arg.lb,
				traceback_arg.rb,
				current_segment_texts,
				suffix_numbers,
				segment_texts
			);
		}
		
		print_segment_texts();
	}
	
	
	void generate_context::calculate_segmentation_lp_copy_segments(
		std::size_t const max_segment_size,
		fseq::segment_text_matrix &segment_texts,
		std::size_t const sequence_count
	)
	{
		struct compare_segments
		{
			std::size_t operator()(fseq::segment_text const &segment) const { return segment.sequence_count(); }
		};
		
		std::vector <std::size_t> count_vector;
		compare_segments cs;
		fseq::segment_text_vector dst_vec, copied_texts;
		
		for (auto &segment_text_vec : segment_texts)
		{
			auto const segment_size(segment_text_vec.size());
			if (segment_size < max_segment_size)
			{
				dst_vec.clear();
				dst_vec.reserve(max_segment_size);
				dst_vec.resize(segment_size);
				copied_texts.clear();
				
				// Sort the current segment by sequence count.
				lb::counting_sort(segment_text_vec, dst_vec, count_vector, cs);
				
				// Copy some of the segments to copied_texts.
				auto const remaining_count(max_segment_size - segment_size);
				for (std::size_t i(0), j(0); j < segment_size; ++j)
				{
					auto const idx(segment_size - j - 1);
					auto segment_text(dst_vec[idx]); // Copy.
					assert(!segment_text.is_copied());
					segment_text.copied_from = idx;
					
					float const proportion(1.0 * segment_text.sequence_count() / sequence_count);
					float const count(proportion * remaining_count);
					auto const copy_count(
						std::min(
							remaining_count - i,
							decltype(remaining_count - i)(std::ceil(count))
						)
					);
					std::fill_n(std::back_inserter(copied_texts), copy_count, segment_text);
					i += copy_count;
					
					if (i == remaining_count)
						break;
				}
				
				// Move the new segments to dst_vec.
				dst_vec.insert(
					dst_vec.end(),
					std::make_move_iterator(copied_texts.begin()),
					std::make_move_iterator(copied_texts.end())
				);
				
				// Store the new value.
				using std::swap;
				swap(segment_text_vec, dst_vec);
			}
		}
		
		print_segment_texts();
	}
	
	
	void generate_context::calculate_segmentation_long_path(std::size_t const lb, std::size_t const rb)
	{
		assert(m_segment_length < rb - lb);
		
		pbwt_context pbwt_ctx(m_sequences, m_alphabet);
		calculate_segmentation_lp_fill_traceback(lb, rb, pbwt_ctx, m_segmentation_traceback);
		
		// TODO: for partitioning, the maximum segment size needs to be calculated over all partitions.
		auto const max_segment_size(m_segmentation_traceback.back().segment_max_size);
		
		calculate_segmentation_lp_create_segments(m_segmentation_traceback, max_segment_size, pbwt_ctx, m_segment_texts);
		calculate_segmentation_lp_copy_segments(max_segment_size, m_segment_texts, pbwt_ctx.size());
	}
	
	
	void generate_context::calculate_segmentation(std::size_t const lb, std::size_t const rb)
	{
		std::cerr << "Calculating the segmentation…" << std::endl;

		if (rb - lb < 2 * m_segment_length)
			calculate_segmentation_short_path(lb, rb);
		else
			calculate_segmentation_long_path(lb, rb);
	}
	
	
	void generate_context::output_segments(lb::file_ostream &stream) const
	{
		// Output format (semi-long form):
		// 1. Segment number
		// 2. Left bound (inclusive)
		// 3. Right bound (exclusive)
		// 4. Segment size
		// 5. Subsequence within the segment
		// 6. List of sequence identifiers in which the subsequence occurs separated by commas.
		// 7. From where the subsequence was copied (for bipartite matching).
		stream << "SEGMENT" "\t" "LB" "\t" "RB" "\t" "SIZE" "\t" "SUBSEQUENCE" "\t" "SEQUENCES" "\t" "COPIED_FROM" "\n";
		for (auto const &tup : boost::combine(m_segmentation_traceback, m_segment_texts))
		{
			auto const &traceback_arg(tup.get <0>());
			auto const &segment_texts(tup.get <1>());
			
			for (auto const &seg_text : segment_texts)
			{
				stream << traceback_arg.lb << '\t' << traceback_arg.rb << '\t' << traceback_arg.segment_size << '\t' << seg_text.text << '\t';
				std::copy(
					seg_text.sequence_indices.cbegin(),
					seg_text.sequence_indices.cend(),
					std::experimental::make_ostream_joiner(stream, ",")
				);
				
				stream
				<< '\t'
				<< (seg_text.is_copied() ? "-" : std::to_string(seg_text.copied_from))
				<< '\n';
			}
		}
		stream << std::flush;
	}
	
	
	void generate_context::match_segments_random_and_output()
	{
		auto const segment_count(m_segment_texts.size());
		auto const segment_size(m_segment_texts.front().size());
		
		for (std::size_t i(0); i < segment_size; ++i)
		{
			for (std::size_t j(0); j < segment_count; ++j)
				std::cout << m_segment_texts[j][i].text;
			std::cout << std::endl;
		}
	}
	
	
	void generate_context::match_segments_sd_and_output()
	{
		// Start tasks for matching the segments.
		auto const segment_count(m_segment_texts.size());
		m_matchings.resize(segment_count - 1);
		for (std::size_t i(1); i < segment_count; ++i)
		{
			m_merge_tasks.emplace_back(new fseq::merge_segments_task(i - 1, m_segment_texts[i - 1], m_segment_texts[i]));
			
			dispatch_group_async(*m_matching_group, *m_queue, ^{
				fseq::matching_vector matchings;
				m_merge_tasks[i - 1]->execute(*m_queue, matchings);
				m_merge_tasks[i - 1].reset();
				
				using std::swap;
				swap(m_matchings[i - 1], matchings);
			});
		}
		
		// Output each founder sequence by following the matched pairs.
		dispatch_group_notify(*m_matching_group, dispatch_get_main_queue(), ^{
			auto const segment_size(m_segment_texts.front().size());
			auto const &leftmost_segment_texts(m_segment_texts[0]);
			
			for (std::size_t i(0); i < segment_size; ++i)
			{
				// Output the leftmost segment texts in order.
				std::size_t lhs(i);
				std::cout << leftmost_segment_texts[lhs].text;
				
				for (std::size_t col(1); col < segment_count; ++col)
				{
					auto const &matchings(m_matchings[col - 1]);
					auto const &segment_texts(m_segment_texts[col]);
					
					// Find the pair and output.
					auto const rhs(matchings[lhs]);
					std::cout << segment_texts[rhs].text;
					
					lhs = rhs;
				}
				
				std::cout << std::endl;
			}
			
			// Exit.
			cleanup();
			exit(EXIT_SUCCESS);
		});
	}
	
	
	void generate_context::output_founders_short_path() const
	{
		if (m_segment_texts.size())
		{
			for (auto const &segment_text : m_segment_texts.front())
				std::cout << segment_text.text << std::endl;
		}
	}
	
	
	void generate_context::prepare()
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
	}
	
	
	void generate_context::load_and_generate(
		char const *input_path,
		fseq::input_format const input_file_format,
		char const *output_segments_path
	)
	{
		lb::file_ostream segments_ostream;
		if (output_segments_path)
			lb::open_file_for_writing(output_segments_path, segments_ostream, false);

		load_input(input_path, input_file_format);
		check_input();
		generate_alphabet();
		
		auto const sequence_length(m_sequences.front().size());
		calculate_segmentation(0, sequence_length);
		
		// Output the segments if needed.
		if (output_segments_path)
			output_segments(segments_ostream);
		
		if (m_segment_texts.size() <= 1)
		{
			output_founders_short_path();
			
			// Exit.
			cleanup();
			exit(EXIT_SUCCESS);
		}
		else
		{
			std::cerr << "Matching the segment contents…" << std::endl;

			switch (m_segment_joining_method)
			{
				case fseq::segment_joining::MATCHING:
					match_segments_sd_and_output();
					break;
				
				case fseq::segment_joining::RANDOM:
					match_segments_random_and_output();
					
					// Exit.
					cleanup();
					exit(EXIT_SUCCESS);
				
				default:
					lb::fail("Unexpected joining method.");
			}
		}
	}
	
	
	void generate_context::print_segment_texts() const
	{
#if PRINT_SEGMENT_TEXTS
		std::size_t vec_idx(0);
		for (auto const &vec : m_segment_texts)
		{
			++vec_idx;
			std::size_t text_idx(0);
			for (auto const &text : vec)
				std::cerr << vec_idx << '.' << ++text_idx << ": " << text << "\n";
		}
		std::cerr << std::flush;
#endif
	}
	
	
	void generate_context::print_segmentation_traceback(
		std::size_t const idx,
		segmentation_dp_arg const &arg
	) const
	{
#if PRINT_SEGMENTATION_TRACEBACK
		std::cerr << "segmentation_traceback_dp[" << idx << "]: " << arg << std::endl;
#endif
	}
}


namespace founder_sequences {
	void generate_founder_sequences(
		char const *input_path,
		fseq::input_format const input_file_format,
		std::size_t const segment_length,
		segment_joining const segment_joining_method,
		char const *output_segments_path,
		bool const use_single_thread
	)
	{
		auto *ctx(new generate_context(segment_length, segment_joining_method, use_single_thread));
	
		ctx->prepare();
		ctx->load_and_generate(input_path, input_file_format, output_segments_path);
	}
}
