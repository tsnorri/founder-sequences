/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef FOUNDER_SEQUENCES_UPDATE_PBWT_TASK_HH
#define FOUNDER_SEQUENCES_UPDATE_PBWT_TASK_HH

#include <founder_sequences/founder_sequences.hh>
#include <vector>


namespace founder_sequences {
	
	class update_pbwt_task final
	{
	protected:
		typedef buffering_pbwt_context::sample_context_type					pbwt_sample_context_type;
		typedef std::vector <std::size_t>									index_vector;
		
	public:
		typedef std::vector <pbwt_sample_type>								pbwt_sample_vector;
		
	protected:
		pbwt_sample_type			m_pbwt_sample;
		pbwt_sample_vector			m_samples;
		index_vector				m_right_bounds;
		std::size_t					m_left_bound{};
		
	public:
		update_pbwt_task() = default;
		
		update_pbwt_task(std::size_t const left_bound, pbwt_sample_type &&pbwt_sample, index_vector &&right_bounds):
			m_pbwt_sample(std::move(pbwt_sample)),
			m_right_bounds(std::move(right_bounds)),
			m_left_bound(left_bound)
		{
		}
		
		pbwt_sample_vector &samples() { return m_samples; }
		pbwt_sample_vector const &samples() const { return m_samples; }
		std::size_t left_bound() const { return m_left_bound; }
		void update_pbwt();
	};
}

#endif
