/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <founder_sequences/segmentation_sp_context.hh>

namespace lb = libbio;


namespace founder_sequences {
	
	void segmentation_sp_context::output_sequence(std::ostream &os, sequence const &seq) const
	{
		auto const bytes(std::as_bytes(seq));
		os.write(reinterpret_cast <char const *>(bytes.data()), bytes.size());
		os << '\n';
	}
	
	
	void segmentation_sp_context::process()
	{
		m_ctx.prepare();
		
		for (std::size_t j(m_lb); j < m_rb; ++j)
		{
			m_ctx.build_prefix_and_divergence_arrays(j);
			m_ctx.update_divergence_value_counts();
			m_ctx.swap_input_and_output();
		}
	}
	
	
	void segmentation_sp_context::output(std::ostream &os, sequence_vector const &sequences) const
	{
		auto const &permutation(m_ctx.input_permutation());
		output_sequence(os, sequences[permutation[0]]);
		
		std::size_t idx(1);
		for (auto const dd : m_ctx.input_divergence() | ranges::view::drop(1))
		{
			if (m_lb < dd)
				output_sequence(os, sequences[permutation[idx]]);
			++idx;
		}
	}
}
