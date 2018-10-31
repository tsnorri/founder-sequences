/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <algorithm>
#include <boost/range/combine.hpp>
#include <iomanip>
#include <iostream>
#include <libbio/assert.hh>
#include <libbio/dispatch_fn.hh>
#include <libbio/file_handling.hh>
#include <libbio/mmap_handle.hh>
#include <range/v3/all.hpp>
#include <regex>

#include "cmdline.h"


namespace ios	= boost::iostreams;
namespace lb	= libbio;


namespace {
	
	typedef ios::stream <ios::array_source>	array_istream;
	typedef std::vector <lb::file_istream>	input_file_vector;
	typedef std::vector <array_istream>		array_istream_vector;
	typedef std::vector <lb::file_ostream>	output_file_vector;
	typedef std::vector <std::string>		file_name_vector;
	
	
	template <typename t_value>
	void open_file_for_writing(
		t_value name,
		lb::file_ostream &stream,
		lb::writing_open_mode const mode
	)
	{
		lb::open_file_for_writing(std::to_string(name), stream, mode);
	}
	
	
	template <>
	void open_file_for_writing(
		std::string const &name,
		lb::file_ostream &stream,
		lb::writing_open_mode const mode
	)
	{
		lb::open_file_for_writing(name, stream, mode);
	}
	
	
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
		for (auto const &tup : ranges::view::zip(input_file_names, input_files))
			lb::open_file_for_reading(std::get <0>(tup), std::get <1>(tup));
	}
	
	
	template <typename t_range>
	void open_output_files(
		t_range const &output_file_names,
		output_file_vector &output_files,
		bool const overwrite
	)
	{
		auto const mode(lb::make_writing_open_mode({
			lb::writing_open_mode::CREATE,
			(overwrite ? lb::writing_open_mode::OVERWRITE : lb::writing_open_mode::NONE)
		}));
		
		for (auto const &tup : ranges::view::zip(output_file_names, output_files))
			open_file_for_writing(std::get <0>(tup), std::get <1>(tup), mode);
	}
	
	
	template <typename t_istream>
	void output_from_streams(
		std::vector <t_istream> &input_streams,
		output_file_vector &output_streams
	)
	{
		lb::for_each(
			boost::combine(input_streams, output_streams),
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
	
	
	template <typename t_istream>
	void process(
		std::vector <t_istream> &input_streams,
		output_file_vector &output_streams,
		lb::file_istream &reference_stream,
		lb::file_istream &identity_column_stream
	)
	{
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
					output_from_streams(input_streams, output_streams);
					break;
			
				// Identity.
				case '1':
					output_from_reference(reference_stream, output_streams, aligned_pos);
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
		
			if (0 == aligned_pos % 1000000)
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
			output_streams,
			[](auto &output) {
				output << std::flush;
			}
		);
	}
	
	
	void process_list_file_and_continue(
		lb::file_istream &reference_stream,
		lb::file_istream &identity_column_stream,
		char const *input_path,
		bool const read_from_cin,
		bool const should_overwrite
	)
	{
		file_name_vector input_file_names;
		file_name_vector output_file_names;
		std::regex output_fname_re(".*?([^/]+)$");
		
		if (read_from_cin)
			read_input_file_names(std::cin, output_fname_re, input_file_names, output_file_names);
		else
		{
			lb::file_istream stream;
			libbio::open_file_for_reading(input_path, stream);
			read_input_file_names(stream, output_fname_re, input_file_names, output_file_names);
		}
		
		// Instantiate streams.
		auto const file_count(input_file_names.size());
		input_file_vector input_files(file_count);
		output_file_vector output_files(file_count);
		open_input_files(input_file_names, input_files);
		open_output_files(output_file_names, output_files, should_overwrite);
		
		process(input_files, output_files, reference_stream, identity_column_stream);
	}
	
	
	void process_text_file_and_continue(
		lb::file_istream &reference_stream,
		lb::file_istream &identity_column_stream,
		char const *input_path,
		bool const read_from_cin,
		bool const should_overwrite
	)
	{
		if (read_from_cin)
			libbio_fail("Memory mapping needed for single-file input.");
		
		lb::mmap_handle handle;
		handle.open(input_path);
		
		// Check the sequence size and that all sequences have the same length.
		auto sv(handle.to_string_view());
		auto const sequence_length(sv.find_first_of('\n'));
		std::size_t sequence_count(0);
		auto const limit(sv.size());
		for (std::size_t i(sequence_length); i < limit; i += (1 + sequence_length))
		{
			libbio_always_assert(sv[i] == '\n');
			++sequence_count;
		}
		
		array_istream_vector inputs(sequence_count);
		std::size_t j(0);
		auto const *data(handle.data());
		for (std::size_t i(sequence_length); i < limit; i+= (1 + sequence_length))
		{
			auto const start(i - sequence_length);
			inputs[j].open(data + start, sequence_length);
			++j;
		}
		
		output_file_vector output_files(sequence_count);
		open_output_files(ranges::view::iota(1, 1 + sequence_count), output_files, should_overwrite);
		
		process(inputs, output_files, reference_stream, identity_column_stream);
	}
}


int main(int argc, char **argv)
{
	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		exit(EXIT_FAILURE);
	
	// Don't sync with stdio.
	std::ios_base::sync_with_stdio(false);
	
	char const *input_path(args_info.input_arg);
	char const *reference_path(args_info.reference_arg);
	char const *identity_column_path(args_info.identity_columns_arg);
	bool const should_overwrite(args_info.overwrite_flag);
	
	bool const read_from_cin('-' == args_info.input_arg[0] && '\0' == args_info.input_arg[1]);
	
	// Open the streams.
	std::cerr << "Opening the files…" << std::endl;
	lb::file_istream reference_stream;
	lb::file_istream identity_column_stream;
	lb::open_file_for_reading(reference_path, reference_stream);
	lb::open_file_for_reading(identity_column_path, identity_column_stream);
	reference_stream.exceptions(reference_stream.exceptions() | std::istream::badbit | std::istream::failbit);
	identity_column_stream.exceptions(identity_column_stream.exceptions() | std::istream::badbit | std::istream::failbit);
	
	switch (args_info.input_format_arg)
	{
		case input_format_arg_text:
			process_text_file_and_continue(reference_stream, identity_column_stream, args_info.input_arg, read_from_cin, should_overwrite);
			break;
		
		case input_format_arg_listMINUS_file:
			process_list_file_and_continue(reference_stream, identity_column_stream, args_info.input_arg, read_from_cin, should_overwrite);
			break;
		
		default:
			libbio_fail("Unexpected input file format");
	}
	
	cmdline_parser_free(&args_info);
	
	return EXIT_SUCCESS;
}
