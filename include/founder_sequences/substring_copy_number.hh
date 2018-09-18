/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef FOUNDER_SEQUENCES_SUBSTRING_COPY_NUMBER_HH
#define FOUNDER_SEQUENCES_SUBSTRING_COPY_NUMBER_HH

#include <cstdint>
#include <ostream>


namespace founder_sequences {

	struct substring_copy_number final
	{
		std::uint32_t	substring_idx{};
		std::uint32_t	copy_number{};
		std::uint32_t	string_idx{};
		
		substring_copy_number() = default;
		
		substring_copy_number(
			std::uint32_t substring_idx_,
			std::uint32_t copy_number_
		):
			substring_idx(substring_idx_),
			copy_number(copy_number_)
		{
		}
	};
	
	typedef std::vector <substring_copy_number>			substring_copy_number_vector;
	typedef std::vector <substring_copy_number_vector>	substring_copy_number_matrix;
	
	
	inline bool operator<(substring_copy_number const &lhs, substring_copy_number const &rhs) { return lhs.copy_number < rhs.copy_number; }
	
	inline std::ostream &operator<<(std::ostream &os, substring_copy_number const &cn) { os << "idx: " << cn.substring_idx << " cn: " << cn.copy_number; return os; }
	
}

#endif
