/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <boost/range/combine.hpp>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <libbio/file_handling.hh>
#include <regex>
#include <set>
#include <vector>

#include "cmdline.h"


namespace lb = libbio;


namespace {
	
	typedef std::vector <lb::file_ostream>			output_file_vector;
	typedef std::vector <std::vector <char>>		buffer_vector;
	typedef std::vector <std::string>				file_name_vector;
	typedef std::vector <std::uint8_t>				identity_column_vector;
	
	
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
	
	
	void create_output_files(file_name_vector const &fnames, bool const overwrite)
	{
		auto const mode(lb::make_writing_open_mode({
			lb::writing_open_mode::CREATE,
			(overwrite ? lb::writing_open_mode::OVERWRITE : lb::writing_open_mode::NONE)
		}));
		
		for (auto const &fname : fnames)
		{
			lb::file_ostream stream;
			lb::open_file_for_writing(fname.c_str(), stream, mode);
		}
	}
	
	
	std::size_t fill_buffers(file_name_vector const &input_file_names, std::size_t const pos, buffer_vector &output_buffers)
	{
		// Fill output buffers from the given files.
		
		std::size_t const count(32 * 1024);
		std::streamsize read_count(-1);
		for (auto const &tup : boost::combine(input_file_names, output_buffers))
		{
			auto const &fname(tup.get <0>());
			auto &buffer(tup.get <1>());
		
			lb::file_istream stream;
			lb::open_file_for_reading(fname.c_str(), stream);
			stream.seekg(pos);
			buffer.resize(count);
		
			if (stream.read(buffer.data(), count))
			{
				read_count = count;
				continue;
			}
		
			// Check if we got a smaller number of characters from the first file.
			if (-1 == read_count)
			{
				read_count = stream.gcount();
				continue;
			}
		
			if (stream.gcount() != read_count)
			{
				std::cerr << "Got an unexpected number of characters from input." << std::endl;
				exit(EXIT_FAILURE);
			}
		}
		
		return read_count;
	}
	
	
	void output_buffer_contents(
		file_name_vector const &output_file_names,
		std::size_t const count,
		buffer_vector const &output_buffers,
		identity_column_vector &skipped_indices
	)
	{
		std::size_t skipped_count(0);
		std::set <char> unique_characters;
		skipped_indices.assign(count, 0);
		for (std::size_t i(0); i < count; ++i)
		{
			unique_characters.clear();
			for (auto const &buffer : output_buffers)
			{
				auto const c(buffer[i]);
				unique_characters.insert(c);
			}
		
			if (unique_characters.size() == 1)
			{
				skipped_indices[i] = true;
				++skipped_count;
			}
		}
	
		if (count == skipped_count)
			return;
	
		for (auto const &tup : boost::combine(output_file_names, output_buffers))
		{
			auto const &fname(tup.get <0>());
			auto const &buffer(tup.get <1>());
		
			lb::file_ostream stream;
			lb::open_file_for_writing(fname.c_str(), stream, lb::writing_open_mode::NONE);
			stream.seekp(0, std::ios_base::end);
		
			for (std::size_t i(0); i < count; ++i)
			{
				if (skipped_indices[i])
					continue;
			
				stream << buffer[i];
			}
		
			stream << std::flush;
		}
	}
}


int main(int argc, char **argv)
{
	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		exit(EXIT_FAILURE);
	
	bool const should_overwrite(args_info.overwrite_flag);
	
	// Don't sync with stdio.
	std::ios_base::sync_with_stdio(false);

	file_name_vector input_file_names;
	file_name_vector output_file_names;
	buffer_vector output_buffers;

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
	
	output_buffers.resize(output_file_names.size());
	
	cmdline_parser_free(&args_info);
	
	create_output_files(output_file_names, should_overwrite);
	
	std::size_t pos(0);
	std::size_t i(0);
	identity_column_vector skipped_indices;
	while (true)
	{
		for (auto &buf : output_buffers)
			buf.clear();

		std::size_t count(fill_buffers(input_file_names, pos, output_buffers));
		if (0 == count)
			break;
		
		pos += count;
		output_buffer_contents(output_file_names, count, output_buffers, skipped_indices);
		
		// Copy the bits as numbers to std::cout.
		std::copy(skipped_indices.cbegin(), skipped_indices.cend(), std::ostream_iterator <int>(std::cout));
		
		++i;
		if (0 == i % 100)
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
			<< "] At position " << pos << "…" << std::endl;
		}
	}
	
	std::cout << std::endl;
	return EXIT_SUCCESS;
}
