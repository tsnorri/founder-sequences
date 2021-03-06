/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef FOUNDER_SEQUENCES_SEGMENTATION_SP_CONTEXT_HH
#define FOUNDER_SEQUENCES_SEGMENTATION_SP_CONTEXT_HH

#include <founder_sequences/segmentation_context.hh>


namespace founder_sequences {
	
	class segmentation_sp_context;
	
	
	struct segmentation_sp_context_delegate : public virtual segmentation_context_delegate
	{
		virtual void context_did_finish_traceback(segmentation_sp_context &ctx) = 0;
	};

	
	// Short path; just count the distinct strings and output.
	class segmentation_sp_context final : public segmentation_context
	{
	protected:
		pbwt_context											m_ctx;
		std::vector <std::pair <std::uint32_t, std::uint32_t>>	m_permutation;
		std::size_t												m_lb{};
		std::size_t												m_rb{};
		std::uint32_t											m_max_segment_size{};
		segmentation_sp_context_delegate						*m_delegate{};

	public:
		segmentation_sp_context(
			segmentation_sp_context_delegate &delegate,
			std::size_t lb,
			std::size_t rb
		):
			m_ctx(delegate.sequences(), delegate.alphabet()),
			m_lb(lb),
			m_rb(rb),
			m_delegate(&delegate)
		{
		}
		
		std::uint32_t max_segment_size() const override { return m_max_segment_size; }
		
		void process();
		void output_founders() const;
		void output_segments() const;
		
	protected:
		void output_sequence(std::ostream &os, sequence const &seq) const;
	};
}

#endif
