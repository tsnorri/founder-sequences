/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef FOUNDER_SEQUENCES_SEGMENTATION_CONTAINER_HH
#define FOUNDER_SEQUENCES_SEGMENTATION_CONTAINER_HH

#include <founder_sequences/founder_sequences.hh>
#include <vector>


namespace founder_sequences {
	
	struct segmentation_container
	{
		std::vector <pbwt_sample_type>	reduced_pbwt_samples;
		segmentation_traceback_vector	reduced_traceback;
		std::uint32_t					max_segment_size{};
	};
}


namespace boost { namespace serialization {
	
	template <typename t_archive>
	void serialize(t_archive &ar, founder_sequences::segmentation_container &container, unsigned int const version)
	{
		ar & container.reduced_pbwt_samples;
		ar & container.reduced_traceback;
		ar & container.max_segment_size;
	}
}}

#endif
