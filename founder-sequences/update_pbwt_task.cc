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
		m_pbwt_sample.set_sample_rate(std::numeric_limits <std::uint64_t>::max());
		for (auto const rb : m_right_bounds)
		{
			m_pbwt_sample.process <lb::pbwt::context_field::DIVERGENCE_VALUE_COUNTS>(rb, [](){});
			
			// Create a sample and copy the fields.
			// FIXME: with pbwt_rmq, the input_permutation will be used as a buffer for RMQ calculation. It does not affect the calculation, though, since at the end of process(), the buffers will be swapped, and the final output permutations will end up as input permutations. API should be added for marking the buffers such that their contents are needed intact after calling process() or something along those lines.
			auto &copied_sample(m_samples.emplace_back());
			copied_sample.set_fields_in_use(
				static_cast <lb::pbwt::context_field>(
					lb::pbwt::context_field::INPUT_PERMUTATION | lb::pbwt::context_field::INPUT_DIVERGENCE
				)
			);
			copied_sample.copy_fields_in_use(m_pbwt_sample);
		}
		
		m_pbwt_sample.set_fields_in_use(lb::pbwt::context_field::NONE);
		m_pbwt_sample.clear_unused_fields();
		m_delegate->task_did_finish(*this);
	}
}
