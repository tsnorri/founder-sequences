/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef FOUNDER_SEQUENCES_GREEDY_MATCHER_HH
#define FOUNDER_SEQUENCES_GREEDY_MATCHER_HH

#include <founder_sequences/founder_sequences.hh>
#include <founder_sequences/matcher.hh>
#include <founder_sequences/substring_copy_number.hh>


namespace founder_sequences { namespace detail {

	struct substring_index_pair final
	{
		std::uint32_t	lhs_idx{};
		std::uint32_t	rhs_idx{};
		std::uint32_t	count{};
		
		substring_index_pair() = default;
		
		substring_index_pair(std::uint32_t lhs_idx_, std::uint32_t rhs_idx_, std::uint32_t count_):
			lhs_idx(lhs_idx_),
			rhs_idx(rhs_idx_),
			count(count_)
		{
		}
		
		substring_index_pair(std::uint32_t lhs_idx_, std::uint32_t rhs_idx_):
			lhs_idx(lhs_idx_),
			rhs_idx(rhs_idx_)
		{
		}
	};
	
	inline bool operator==(substring_index_pair const &lhs, substring_index_pair const &rhs)
	{
		return lhs.lhs_idx == rhs.rhs_idx && lhs.rhs_idx == rhs.rhs_idx;
	}
	
	inline std::ostream &operator<<(std::ostream &os, substring_index_pair const &pair)
	{
		os << '(' << pair.lhs_idx << ", " << pair.rhs_idx << "; " << pair.count << ')';
		return os;
	}
}}


namespace founder_sequences
{
	class greedy_matcher;
	
	
	struct greedy_matcher_delegate : public virtual matcher_delegate
	{
		virtual void matcher_did_finish(greedy_matcher &matcher) = 0;
	};
	
	
	class greedy_matcher final : public matcher
	{
	protected:
		typedef std::tuple <substring_copy_number_vector, pbwt_sample_type, permutation_vector>	matching_tuple;
		
	protected:
		greedy_matcher_delegate					*m_delegate{};
		substring_copy_number_matrix const		*m_substrings_to_output{};
		std::size_t								m_permutation_max{};
		std::size_t								m_permutation_bits_needed{};
		
	public:
		greedy_matcher(
			greedy_matcher_delegate &delegate,
			substring_copy_number_matrix const &substrings_to_output,
			std::size_t const permutation_max,
			std::size_t const permutation_bits_needed
		):
			m_delegate(&delegate),
			m_substrings_to_output(&substrings_to_output),
			m_permutation_max(permutation_max),
			m_permutation_bits_needed(permutation_bits_needed)
		{
		}
		
		void match() override;
		void output_segments(std::ostream &stream, sequence_vector const &sequences) override;
		
	protected:
		auto create_index_pairs(
			matching_tuple const &lhs,
			matching_tuple const &rhs,
			sdsl::int_vector <0> const &rhs_matching,
			std::vector <detail::substring_index_pair> &index_pairs,
			std::vector <detail::substring_index_pair> &index_pairs_buffer,
			std::vector <std::uint32_t> &lhs_unused_substring_numbers,
			sdsl::int_vector <0> &to_lhs_substring,
			sdsl::int_vector <0> &to_rhs_string
		) const -> std::pair <std::size_t, std::size_t>;
		
		void create_matching(
			std::vector <detail::substring_index_pair> const &index_pairs,
			std::uint64_t const matching_max,
			sdsl::int_vector <0> &lhs_matching,
			sdsl::int_vector <0> &rhs_matching
		) const;
	};
}

#endif
