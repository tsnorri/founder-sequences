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
		m_ctx.process <lb::pbwt::context_field::DIVERGENCE_VALUE_COUNTS>(m_rb, [](std::size_t const){});
		
		m_max_segment_size = m_ctx.unique_substring_count_lhs(m_lb);
		m_ctx.unique_substring_count_idxs_lhs(0, m_permutation);
		m_delegate->context_did_finish_traceback(*this);
	}
	
	
	void segmentation_sp_context::output_founders() const
	{
		auto const &sequences(m_delegate->sequences());
		auto &stream(m_delegate->sequence_output_stream());
		
		for (auto const &pair : m_permutation)
			output_sequence(stream, sequences[pair.first]);
	}
	
	
	void segmentation_sp_context::output_segments() const
	{
		auto &stream(m_delegate->segments_output_stream());
		stream << "SEQUENCE" "\n";
		for (auto const &pair : m_permutation)
			stream << pair.second << '\n';
	}
}
