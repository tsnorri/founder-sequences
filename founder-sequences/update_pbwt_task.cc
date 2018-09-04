/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <founder_sequences/update_pbwt_task.hh>


namespace founder_sequences {

	void update_pbwt_task::update_pbwt()
	{
		// Take the next right bound, update the sample up to it.
		for (auto const rb : m_right_bounds)
		{
			++m_pbwt_sample.index;
			while (m_pbwt_sample.index < rb)
			{
				m_pbwt_sample.context.build_prefix_and_divergence_arrays(m_pbwt_sample.index);
				m_pbwt_sample.context.update_divergence_value_counts();
				m_pbwt_sample.context.swap_input_and_output();
				++m_pbwt_sample.index;
			}
			
			// Copy the updated sample.
			m_samples.emplace_back(m_pbwt_sample);
		}
		
		m_pbwt_sample.context.clear();
	}
}
