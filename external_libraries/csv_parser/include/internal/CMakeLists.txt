add_library(csv_parser STATIC "")

target_sources(csv_parser
	PRIVATE
		compatibility.hpp
		constants.hpp
		csv_format.hpp
		csv_format.cpp
		csv_reader.hpp
		csv_reader.cpp
		csv_reader_internals.hpp
        csv_reader_internals.cpp
		csv_reader_iterator.cpp
		csv_row.hpp
		csv_row.cpp
		csv_row_json.cpp
		csv_stat.cpp
		csv_stat.hpp
		csv_utility.cpp
		csv_utility.hpp
		csv_writer.hpp
		data_type.h
		row_buffer.hpp
		row_buffer.cpp
)

set_target_properties(csv_parser PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(csv_parser PRIVATE Threads::Threads)

set(CSV_PARSER_LIBRARIES csv_parser)

message(STATUS "csv_parser : ${CSV_INCLUDE_DIR}")
