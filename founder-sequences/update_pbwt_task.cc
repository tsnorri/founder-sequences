/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <founder_sequences/update_pbwt_task.hh>

namespace lb = libbio;


namespace founder_sequences {

	void update_pbwt_task::update_pbwt()
	{
		// Take the next right bound, update the sample up to it.
		for (auto const rb : m_right_bounds)
		{
			while (m_pbwt_sample.rb < rb)
			{
				m_pbwt_sample.context.build_prefix_and_divergence_arrays(m_pbwt_sample.rb);
				m_pbwt_sample.context.update_divergence_value_counts();
				m_pbwt_sample.context.swap_input_and_output();
				++m_pbwt_sample.rb;
			}
			
			// Copy the updated sample.
			auto &copied_sample(m_samples.emplace_back(m_pbwt_sample));
			
			// Release memory.
			copied_sample.context.clear(
				static_cast <lb::pbwt::pbwt_context_field>(
					lb::pbwt::pbwt_context_field::ALL & ~lb::pbwt::pbwt_context_field::INPUT_DIVERGENCE
				)
			);
		}
		
		m_pbwt_sample.context.clear();
	}
}
