/*
 * Copyright (c) 2018 Dmitry Kosolobov, Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef FOUNDER_SEQUENCES_RMQ_HH
#define FOUNDER_SEQUENCES_RMQ_HH

#include <algorithm>
#include <cstdint>
#include <sdsl/bits.hpp>
#include <vector>


namespace founder_sequences {
	
	template <
		typename t_values,
		typename t_cmp = std::less <typename t_values::value_type>,
		unsigned int t_block_size = 32
	>
	class rmq {
		static_assert((t_block_size & (t_block_size - 1)) == 0, "Block size must be a power of two");
	
	protected:
		std::vector <
			std::vector <
				typename t_values::size_type
			>
		>																	m_precalc;
		t_values const														*m_values{nullptr};
		t_cmp																m_cmp;
	
	public:
		rmq() = default;
		
		rmq(t_values const &values, t_cmp cmp):
			m_precalc(1),
			m_values(&values),
			m_cmp(cmp)
		{
		}
		
		explicit rmq(t_values const &values):
			m_precalc(1),
			m_values(&values)
		{
		}
		
		void set_values(t_values const &values) { m_values = &values; }
		void update(std::size_t last_idx);
		std::size_t operator()(std::size_t beg, std::size_t end) const;
	
	protected:
		std::size_t naive_min(std::size_t first, std::size_t last) const;
	};
	
	
	// Update the RMQ data structure after appending a value to m_values.
	template <typename t_values, typename t_cmp, unsigned int t_block_size>
	void rmq <t_values, t_cmp, t_block_size>::update(std::size_t const last_idx)
	{
		if (((1 + last_idx) & (t_block_size - 1)) != 0)
			return;
		
		m_precalc.resize(1 + (1 + last_idx) / t_block_size);
		std::size_t bnum = 1 + last_idx / t_block_size;
		std::size_t new_smp = naive_min((bnum - 1) * t_block_size, bnum * t_block_size);
		auto const &new_val = (*m_values)[new_smp];
		assert(0 < m_precalc.size());
		m_precalc[0].push_back(new_smp);
		for (std::size_t pow2 = 1; (1u << pow2) <= bnum; ++pow2)
		{
			assert(pow2 < m_precalc.size());
			assert(bnum - (1u << pow2) < m_precalc[pow2 - 1].size());
			std::size_t smp1 = m_precalc[pow2 - 1][bnum - (1u << pow2)];
			std::size_t smp2 = m_precalc[pow2 - 1][bnum - ((1u << pow2) - 1) - 1];
			std::size_t smp = (m_cmp((*m_values)[smp2], (*m_values)[smp1]) ? smp2 : smp1);
			m_precalc[pow2].push_back(m_cmp(new_val, (*m_values)[smp]) ? new_smp : smp);
		}
	}

	// Find the position of a minimal element in the range [beg..end).
	template <typename t_values, typename t_cmp, unsigned int t_block_size>
	std::size_t rmq <t_values, t_cmp, t_block_size>::operator()(std::size_t beg, std::size_t end) const
	{
		std::size_t beg_block = beg / t_block_size + 1;
		std::size_t end_block = end / t_block_size;
	
		if (beg_block >= end_block)
			return naive_min(beg, end);
	
		std::size_t const pow2 = sdsl::bits::hi(end_block - beg_block);
		std::size_t const smp1 = m_precalc[pow2][beg_block];
		std::size_t const smp2 = m_precalc[pow2][end_block - (1u << pow2)];
		std::size_t smp = (m_cmp((*m_values)[smp2], (*m_values)[smp1]) ? smp2 : smp1);
		std::size_t const left_smp = naive_min(beg, beg_block * t_block_size);
		smp = (m_cmp((*m_values)[left_smp], (*m_values)[smp]) ? left_smp : smp);
	
		if (end == end_block * t_block_size)
			return smp;
	
		std::size_t right_smp = naive_min(end_block * t_block_size, end);
		return (m_cmp((*m_values)[right_smp], (*m_values)[smp]) ? right_smp : smp);
	}
	
	
	template <typename t_values, typename t_cmp, unsigned int t_block_size>
	std::size_t rmq <t_values, t_cmp, t_block_size>::naive_min(std::size_t first, std::size_t last) const
	{
		assert(nullptr != m_values);
		assert(first < m_values->size());
		assert(last <= m_values->size());
		
		auto const begin(m_values->cbegin());
		auto const it(std::min_element(begin + first, begin + last, m_cmp));
		return std::distance(begin, it);
	}
}

#endif
