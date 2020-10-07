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
#include "utilities/string_utilities.h"
#include "includes/kratos_filesystem.h"
#include "csv-parser/include/csv.hpp"

namespace Kratos {

namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(CSVRowInterator, KratosExternalLibrariesFastSuite)
{
    std::string csv_string = "A,B,C\r\n" // Header row
        "123,234,345\r\n"
        "1,2,3\r\n"
        "1,2,3";
    
    csv::CSVRow row;
    auto rows = csv::parse(csv_string);
    rows.read_row(row);

    // Forwards and Backwards Iterators
    {
        // Forwards
        KRATOS_CHECK_EQUAL(row.begin()->get<int>(), 123);
        KRATOS_CHECK_STRING_EQUAL((row.end() - 1)->get<>(), "345");

        std::size_t i = 0;
        for (auto it = row.begin(); it != row.end(); ++it) {
            if (i == 0) {
                KRATOS_CHECK_EQUAL(it->get<>(), "123");
            } else if (i == 1) {
                KRATOS_CHECK_EQUAL(it->get<>(), "234");
            } else {
                KRATOS_CHECK_EQUAL(it->get<>(), "345");
            }

            i++;
        }

        // Backwards
        KRATOS_CHECK_EQUAL(row.rbegin()->get<int>(), 345);
        KRATOS_CHECK_STRING_EQUAL((row.rend() - 1)->get<>(), "123");
    }

    // Iterator Arithmetic
    {
        KRATOS_CHECK_EQUAL(row.begin()->get<int>(), 123);
        KRATOS_CHECK_STRING_EQUAL((row.end() - 1)->get<>(), "345");

        auto row_start = row.begin();
        KRATOS_CHECK_EQUAL(*(row_start + 1), "234");
        KRATOS_CHECK_EQUAL(*(row_start + 2), "345");

    }

    // Post-Increment Iterator
    {
        auto it = row.begin();

        KRATOS_CHECK_EQUAL(it++->get<int>(), 123);
        KRATOS_CHECK_EQUAL(it->get<int>(), 234);

        KRATOS_CHECK_EQUAL(it--->get<int>(), 234);
        KRATOS_CHECK_EQUAL(it->get<int>(), 123);
    }

    // Range Based For
    {
        std::size_t i = 0;
        for (auto& field : row) {
            if (i == 0) {
                KRATOS_CHECK_EQUAL(field.get<>(), "123");
            } else if (i == 1) {
                KRATOS_CHECK_EQUAL(field.get<>(), "234");
            } else {
                KRATOS_CHECK_EQUAL(field.get<>(), "345");
            }

            i++;
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(BasicCSVReaderIteratorTest, KratosExternalLibrariesFastSuite)
{
    // A file with 100 rows and columns A, B, ... J
    // where every value in the ith row is the number i
    const std::string working_dir = StringUtilities::ErasePartialString(__FILE__, "test_csv_iterator.cpp");
    csv::CSVReader reader(Kratos::FilesystemExtensions::JoinPaths({working_dir, "data/fake_data/ints.csv"}));
    std::vector<std::string> col_names = {
        "A", "B", "C", "D", "E", "F", "G", "H", "I", "J"
    };
    std::size_t i = 1;

    for (auto it = reader.begin(); it != reader.end(); ++it) {
        KRATOS_CHECK_EQUAL((*it)[0].get<std::size_t>(), i);
        i++;
    }
}

KRATOS_TEST_CASE_IN_SUITE(CSVReaderIteratorstdmax_elem, KratosExternalLibrariesFastSuite)
{
    // The first is such that each value in the ith row is the number i
    // There are 100 rows
    // The second file is a database of California state employee salaries
    const std::string working_dir = StringUtilities::ErasePartialString(__FILE__, "test_csv_iterator.cpp");
    csv::CSVReader r1(Kratos::FilesystemExtensions::JoinPaths({working_dir, "data/fake_data/ints.csv"}));
    csv::CSVReader r2(Kratos::FilesystemExtensions::JoinPaths({working_dir, "data/real_data/2015_StateDepartment.csv"}));

    // Find largest number
    auto int_finder = [](csv::CSVRow& left, csv::CSVRow& right) {
        return (left["A"].get<int>() < right["A"].get<int>());
    };

    auto max_int = std::max_element(r1.begin(), r2.end(), int_finder);

    // Find highest salary
    auto wage_finder = [](csv::CSVRow& left, csv::CSVRow& right) {
        return (left["Total Wages"].get<double>() < right["Total Wages"].get<double>());
    };

    auto max_wage = std::max_element(r2.begin(), r2.end(), wage_finder);

    KRATOS_CHECK_EQUAL((*max_int)["A"], 100);
    KRATOS_CHECK_EQUAL((*max_wage)["Total Wages"], "812064.87");
}

} // namespace Testing.
} // namespace Kratos.
