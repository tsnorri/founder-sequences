/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <founder_sequences/update_pbwt_task.hh>

namespace lb = libbio;


namespace founder_sequences {

	void update_pbwt_task::execute()
	{
		// Take the next right bound, update the sample up to it.
		for (auto const rb : m_right_bounds)
		{
			m_pbwt_sample.context.process <lb::pbwt::context_field::DIVERGENCE_VALUE_COUNTS>(rb, [](auto const idx){});
			
			// Create a sample and copy the fields.
			auto &copied_sample(m_samples.emplace_back(rb));
			copied_sample.context.set_fields_in_use(
				static_cast <lb::pbwt::context_field>(
					lb::pbwt::context_field::INPUT_PERMUTATION | lb::pbwt::context_field::INPUT_DIVERGENCE
				)
			);
			copied_sample.context.copy_fields_in_use(m_pbwt_sample.context);
		}
		
		m_pbwt_sample.context.set_fields_in_use(lb::pbwt::context_field::NONE);
		m_pbwt_sample.context.clear_unused_fields();
		m_delegate->task_did_finish(*this);
	}
}
