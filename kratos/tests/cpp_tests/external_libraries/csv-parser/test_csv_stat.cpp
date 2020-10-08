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

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/kratos_filesystem.h"
#include "utilities/string_utilities.h"
#include "csv-parser/csv.hpp"

namespace Kratos {

namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(CalculatingStatisticsfromDirectInput, KratosExternalLibrariesFastSuite )
{
    std::string int_str;
    std::string int_list = "";
    for (int i = 1; i < 101; i++) {
        int_str = std::to_string(i);
        int_list += int_str + "," + int_str + "," + int_str + "\r\n";
    }

    // Expected results
    csv::CSVFormat format;
    format.column_names({ "A", "B", "C" });

    csv::CSVStat reader(format);
    reader.feed(int_list);
    reader.end_feed();

    std::vector<long double> means = { 50.5, 50.5, 50.5 };
    std::vector<long double> mins = { 1, 1, 1 };
    std::vector<long double> maxes = { 100, 100, 100 };

    KRATOS_CHECK_EQUAL( reader.get_mins(), mins );
    KRATOS_CHECK_EQUAL( reader.get_maxes(), maxes );
    KRATOS_CHECK_EQUAL( reader.get_mean(), means );
    KRATOS_CHECK_EQUAL( ceil(reader.get_variance()[0]), 842 );

    // Make sure all integers between 1 and 100 have a count of 1
    for (int i = 1; i < 101; i++)
        KRATOS_CHECK_EQUAL( reader.get_counts()[0][std::to_string(i)], 1 );

    // Confirm column at pos 0 has 100 integers (type 2)
    KRATOS_CHECK_EQUAL( reader.get_dtypes()[0][csv::DataType::CSV_INT8], 100 );
}
} // namespace Testing.
} // namespace Kratos.
