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
#include "csv-parser/include/csv.hpp"

namespace Kratos {

namespace Testing {

// Construct a CSVRow and assert that its interface works as expected
KRATOS_TEST_CASE_IN_SUITE(CSVRowTest, KratosExternalLibrariesFastSuite) {
    std::string csv_string = "A,B,C,D\r\n"
                  "Col1,Col2,Col3,Col4";

    csv::CSVRow row;
    auto rows = csv::parse(csv_string);
    rows.read_row(row);

    bool error_caught = false;

    // size() Check
    {
        KRATOS_CHECK_EQUAL(row.size(), 4);
    }

    // operator[]
    {
        KRATOS_CHECK_EQUAL(row[1], "Col2");
        KRATOS_CHECK_EQUAL(row["B"], "Col2");

        KRATOS_CHECK_EQUAL(row[2], "Col3");
        KRATOS_CHECK_EQUAL(row["C"], "Col3");
    }

    // operator[] Out of Bounds
    {
        try {
            row[4].get<>();
        }
        catch (std::runtime_error&) {
            error_caught = true;
        }

        KRATOS_CHECK(error_caught);
    }

    // operator[] Access Non-Existent Column
    {
        try {
            row["Col5"].get<>();
        }
        catch (std::runtime_error&) {
            error_caught = true;
        }

        KRATOS_CHECK(error_caught);
    }

    // Content Check
    {
        KRATOS_CHECK_EQUAL(std::vector<std::string>(row),
            std::vector<std::string>({ "Col1", "Col2", "Col3", "Col4" }));
    }

    /** Allow get_sv() to be used with a const CSVField
     *
     *  See: https://github.com/vincentlaucsb/csv-parser/issues/86
     *
     */
    // get_sv() Check
    {
        std::vector<std::string> content;

        for (const auto& field : row) {
            content.push_back(std::string(field.get_sv()));
        }

        KRATOS_CHECK_EQUAL(std::vector<std::string>(row),
            std::vector<std::string>({ "Col1", "Col2", "Col3", "Col4" }));
    }
}

// Integration test for CSVRow/CSVField
KRATOS_TEST_CASE_IN_SUITE(CSVFieldoperatorEqual, KratosExternalLibrariesFastSuite) {
    std::string csv_string = "A,B,C,D\r\n"
                  "1,2,3,3.14";

    csv::CSVRow row;
    auto rows = csv::parse(csv_string);
    rows.read_row(row);

    KRATOS_CHECK_EQUAL(row["A"], 1);
    KRATOS_CHECK_EQUAL(row["B"], 2);
    KRATOS_CHECK_EQUAL(row["C"], 3);
    KRATOS_CHECK(csv::internals::is_equal(row["D"].get<long double>(), 3.14L));
}

} // namespace Testing.
} // namespace Kratos.
