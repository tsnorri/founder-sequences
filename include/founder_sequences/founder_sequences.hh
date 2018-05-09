/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef FOUNDER_SEQUENCES_FOUNDER_SEQUENCES_HH
#define FOUNDER_SEQUENCES_FOUNDER_SEQUENCES_HH

#include <memory>
#include <string>
#include <vector>


namespace founder_sequences {
	
	struct segment_text
	{
		std::string					text;
		std::vector <std::size_t>	sequence_indices;
		std::size_t					copied_from{SIZE_MAX};
		
		bool is_copied() const { return SIZE_MAX != copied_from; }
		std::size_t first_sequence_index() const { return sequence_indices.front(); }
		std::size_t sequence_count() const { return sequence_indices.size(); }
		std::size_t row_number(std::size_t const row) const { return (SIZE_MAX == copied_from ? row : copied_from); }
	};
	
	std::ostream &operator<<(std::ostream &os, segment_text const &seg_text);
	
	typedef std::vector <segment_text>			segment_text_vector;
	typedef std::vector <segment_text_vector>	segment_text_matrix;
	
	
	typedef std::vector <std::uint32_t>			matching_vector;
	
	
	enum class segment_joining : uint8_t {
		MATCHING	= 0,
		RANDOM
	};
	
	enum class input_format : uint8_t {
		FASTA		= 0,
		LIST_FILE
	};
	
	
	typedef std::vector <std::vector <std::uint8_t>>	sequence_vector;
	
	template <typename t_vector_source>
	class reader_cb
	{
	protected:
		sequence_vector	*m_sequences{nullptr};
		
	public:
		reader_cb(sequence_vector &vec):
			m_sequences(&vec)
		{
		}
		
		void handle_sequence(
			std::unique_ptr <typename t_vector_source::vector_type> &seq_ptr,
			std::size_t const &seq_length,
			t_vector_source &vector_source
		)
		{
			auto &seq(m_sequences->emplace_back(*seq_ptr));
			seq.resize(seq_length);
			seq_ptr.reset();
			vector_source.set_vector_length(seq_length);
			// Don't return the vector to the source.
		}
		
		void start() {}
		void finish() {}
	};
	
	
	template <typename t_vector_source>
	class fasta_reader_cb final : public reader_cb <t_vector_source>
	{
	public:
		using reader_cb <t_vector_source>::reader_cb;
		
		void handle_sequence(
			std::string const &identifier,
			std::unique_ptr <typename t_vector_source::vector_type> &seq_ptr,
			std::size_t const &seq_length,
			t_vector_source &vector_source
		)
		{
			reader_cb <t_vector_source>::handle_sequence(seq_ptr, seq_length, vector_source);
		}
	};
	
	
	template <typename t_vector_source>
	class line_reader_cb final : public reader_cb <t_vector_source>
	{
	public:
		using reader_cb <t_vector_source>::reader_cb;
		
		void handle_sequence(
			uint32_t line,
			std::unique_ptr <typename t_vector_source::vector_type> &seq_ptr,
			std::size_t const &seq_length,
			t_vector_source &vector_source
		)
		{
			reader_cb <t_vector_source>::handle_sequence(seq_ptr, seq_length, vector_source);
		}
	};
}

#endif
