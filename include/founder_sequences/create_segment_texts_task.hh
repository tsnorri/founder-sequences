/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef FOUNDER_SEQUENCES_CREATE_SEGMENT_TEXTS_TASK_HH
#define FOUNDER_SEQUENCES_CREATE_SEGMENT_TEXTS_TASK_HH

#include <founder_sequences/founder_sequences.hh>
#include <founder_sequences/segment_text.hh>
#include <founder_sequences/substring_copy_number.hh>


namespace founder_sequences {
	
	class create_segment_texts_task;
	
	
	class create_segment_texts_task final : public task
	{
	protected:
		pbwt_sample_type const						*m_pbwt_sample{};
		substring_copy_number_vector const			*m_substring_copy_numbers{};
		segment_text_vector							*m_segment_texts{};
		std::size_t									m_max_segment_size{};
		std::size_t									m_seq_count{};
		
	public:
		create_segment_texts_task() = default;
		
		
		create_segment_texts_task(
			pbwt_sample_type const &sample,
			substring_copy_number_vector const &substring_copy_numbers,
			std::size_t const max_segment_size,
			std::size_t const seq_count,
			segment_text_vector &segment_texts
		):
			m_pbwt_sample(&sample),
			m_substring_copy_numbers(&substring_copy_numbers),
			m_segment_texts(&segment_texts),
			m_max_segment_size(max_segment_size),
			m_seq_count(seq_count)
		{
		}
		
		void execute() override;
	};
}

#endif
