//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vincent La (This an adaptation of https://github.com/vincentlaucsb/csv-parser)
//                   Vicente Mataix Ferrandiz
//
//

// System includes
#include <stdio.h> // remove()
#include <sstream>

// External includes

// Project includes
#include "testing/testing.h"
#include "utilities/string_utilities.h"
#include "includes/kratos_filesystem.h"
#include "csv-parser/csv.hpp"

namespace Kratos {

namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(PreventColumnNamesFromBeingOverwritten, KratosExternalLibrariesFastSuite)
{
    const std::string working_dir = StringUtilities::ErasePartialString(__FILE__, "test_read_csv_file.cpp");

    std::vector<std::string> column_names = { "A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10" };

    // Test against a variety of different csv::CSVFormat objects
    std::vector<csv::CSVFormat> formats = {};
    formats.push_back(csv::CSVFormat::guess_csv());
    formats.push_back(csv::CSVFormat());
    formats.back().delimiter(std::vector<char>({ ',', '\t', '|'}));
    formats.push_back(csv::CSVFormat());
    formats.back().delimiter(std::vector<char>({ ',', '~'}));

    for (auto& format_in : formats) {
        // Set up the csv::CSVReader
        format_in.column_names(column_names);
        csv::CSVReader reader(Kratos::FilesystemExtensions::JoinPaths({working_dir, "data/fake_data/ints_comments.csv"}), format_in);

        // Assert that column names weren't overwritten
        csv::CSVFormat format_out = reader.get_format();
        KRATOS_CHECK_EQUAL(reader.get_col_names(), column_names);
        KRATOS_CHECK_EQUAL(format_out.get_delim(), ',');
        KRATOS_CHECK_EQUAL(format_out.get_header(), 5);
    }
}

KRATOS_TEST_CASE_IN_SUITE(NonExistentCSV, KratosExternalLibrariesFastSuite)
{
    // Make sure attempting to parse a non-existent CSV throws an error
    bool error_caught = false;

    try {
        csv::CSVReader reader("lochness.csv");
    }
    catch (std::runtime_error& err) {
        error_caught = true;
        KRATOS_CHECK_EQUAL(err.what(), std::string("Cannot open file lochness.csv"));
    }

    KRATOS_CHECK(error_caught);
}

KRATOS_TEST_CASE_IN_SUITE(read_rowCSVField_Easy, KratosExternalLibrariesFastSuite)
{
    // Test that integers are type-casted properly
    const std::string working_dir = StringUtilities::ErasePartialString(__FILE__, "test_read_csv_file.cpp");
    csv::CSVReader reader(Kratos::FilesystemExtensions::JoinPaths({working_dir, "data/fake_data/ints.csv"}));
    csv::CSVRow row;

    while (reader.read_row(row)) {
        for (std::size_t i = 0; i < row.size(); i++) {
            KRATOS_CHECK(row[i].is_int());
            KRATOS_CHECK_LESS_EQUAL(row[i].get<int>(), 100);
        }
    }
}

} // namespace Testing.
} // namespace Kratos.
