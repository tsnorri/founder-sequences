/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <algorithm>
#include <boost/range/combine.hpp>
#include <iomanip>
#include <iostream>
#include <libbio/dispatch_fn.hh>
#include <libbio/file_handling.hh>
#include <range/v3/all.hpp>
#include <regex>

#include "cmdline.h"


namespace lb = libbio;


namespace {
	
	typedef std::vector <lb::file_istream>	input_file_vector;
	typedef std::vector <lb::file_ostream>	output_file_vector;
	typedef std::vector <std::string>		file_name_vector;
	
	
	void read_input_file_names(
		std::istream &istream,
		std::regex const &output_fname_re,
		file_name_vector &input_file_names,
		file_name_vector &output_file_names
	)
	{
		std::string fname;
		while (std::getline(istream, fname))
		{
			input_file_names.push_back(fname);
			
			std::smatch fname_match;
			std::regex_search(fname, fname_match, output_fname_re);
			auto out_fname(fname_match[1].str());
			output_file_names.emplace_back(std::move(out_fname));
		}
	}
	
	
	void open_input_files(
		file_name_vector const &input_file_names,
		input_file_vector &input_files
	)
	{
		for (auto const &tup : boost::combine(input_file_names, input_files))
			lb::open_file_for_reading(tup.get <0>().c_str(), tup.get <1>());
	}
	
	
	void open_output_files(
		file_name_vector const &output_file_names,
		output_file_vector &output_files,
		bool const overwrite
	)
	{
		auto const mode(lb::make_writing_open_mode({
			lb::writing_open_mode::CREATE,
			(overwrite ? lb::writing_open_mode::OVERWRITE : lb::writing_open_mode::NONE)
		}));
		
		for (auto const &tup : boost::combine(output_file_names, output_files))
			lb::open_file_for_writing(tup.get <0>().c_str(), tup.get <1>(), mode);
	}
	
	
	void output_from_streams(
		input_file_vector &input_files,
		output_file_vector &output_files
	)
	{
		auto combined(boost::combine(input_files, output_files));
		lb::for_each(
			boost::combine(input_files, output_files),
			[](auto const &tup) {
				char c{};
				auto &input(tup.template get <0>());
				auto &output(tup.template get <1>());
				input >> std::noskipws >> c;
				output << c;
			}
		);
	}
	
	
	void output_from_reference(
		std::istream &reference,
		output_file_vector &output_files,
		std::size_t const aligned_pos
	)
	{
		char c{};
		reference.seekg(aligned_pos);
		reference >> std::noskipws >> c;
		lb::for_each(
			output_files,
			[c](auto &output) {
				output << c;
			}
		);
	}
}


int main(int argc, char **argv)
{
	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		exit(EXIT_FAILURE);
	
	// Don't sync with stdio.
	std::ios_base::sync_with_stdio(false);
	
	char const *input_list_path(args_info.input_arg);
	char const *reference_path(args_info.reference_arg);
	char const *identity_column_path(args_info.identity_columns_arg);
	bool const should_overwrite(args_info.overwrite_flag);
	
	file_name_vector input_file_names;
	file_name_vector output_file_names;
	
	// Read the input file names and process them.
	std::regex output_fname_re(".*?([^/]+)$");
	if ('-' == args_info.input_arg[0] && '\0' == args_info.input_arg[1])
		read_input_file_names(std::cin, output_fname_re, input_file_names, output_file_names);
	else
	{
		lb::file_istream stream;
		libbio::open_file_for_reading(args_info.input_arg, stream);
		read_input_file_names(stream, output_fname_re, input_file_names, output_file_names);
	}
	
	// Instantiate streams.
	auto const file_count(input_file_names.size());
	input_file_vector input_files(file_count);
	output_file_vector output_files(file_count);
	lb::file_istream reference_stream;
	lb::file_istream identity_column_stream;
	
	// Open the streams.
	std::cerr << "Opening the files…" << std::endl;
	open_input_files(input_file_names, input_files);
	open_output_files(output_file_names, output_files, should_overwrite);
	lb::open_file_for_reading(reference_path, reference_stream);
	lb::open_file_for_reading(identity_column_path, identity_column_stream);
	
	cmdline_parser_free(&args_info);
	
	// Read from the inputs and combine.
	std::cerr << "Handling the input…" << std::endl;
	std::size_t aligned_pos(0);
	char is_identity{};
	while (identity_column_stream >> std::noskipws >> is_identity)
	{
		switch (is_identity)
		{
			// Not identity.
			case '0':
				output_from_streams(input_files, output_files);
				break;
			
			// Identity.
			case '1':
				output_from_reference(reference_stream, output_files, aligned_pos);
				break;
			
			// End of file.
			case '\n':
				goto exit_loop;
			
			// Unexpected character.
			default:
				throw std::runtime_error("Unexpected character");
				break;
		}
		
		++aligned_pos;
		
		if (0 == aligned_pos % 10000)
		{
			auto const time(std::chrono::system_clock::now());
			auto const ct(std::chrono::system_clock::to_time_t(time));
			struct tm ctm;
			localtime_r(&ct, &ctm);
			std::cerr
			<< '['
			<< std::setw(2) << std::setfill('0') << ctm.tm_hour << ':'
			<< std::setw(2) << std::setfill('0') << ctm.tm_min << ':'
			<< std::setw(2) << std::setfill('0') << ctm.tm_sec
			<< "] At position " << aligned_pos << "…" << std::endl;
		}
	}
	
exit_loop:
	lb::for_each(
		output_files,
		[](auto &output) {
			output << std::flush;
		}
	);
	
	return EXIT_SUCCESS;
}
