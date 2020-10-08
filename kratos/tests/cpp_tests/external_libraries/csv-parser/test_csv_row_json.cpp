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
#include "csv-parser/csv.hpp"

namespace Kratos {

namespace Testing {

/** Construct a CSVRow object for testing given column names and CSV fields */
csv::CSVRow make_csv_row(std::vector<std::string> data, std::vector<std::string> col_names) {
    // Concatenate vector or strings into one large string

    std::stringstream raw_csv;
    auto writer = csv::make_csv_writer(raw_csv);
    writer << col_names;
    writer << data;

    csv::CSVReader reader;
    reader.feed(raw_csv.str());
    reader.end_feed();

    csv::CSVRow row;
    reader.read_row(row);

    return row;
}

KRATOS_TEST_CASE_IN_SUITE(json_escape_stringTest, KratosExternalLibrariesFastSuite)
{
    using csv::internals::json_escape_string;

    // Assert that special characters are escaped properly
    KRATOS_CHECK_EQUAL(json_escape_string("Quote\"Quote"), "Quote\\\"Quote");
    KRATOS_CHECK_EQUAL(json_escape_string("RSolidus\\RSolidus"), "RSolidus\\\\RSolidus");
    KRATOS_CHECK_EQUAL(json_escape_string("Backspace\bBackspace"), "Backspace\\bBackspace");
    KRATOS_CHECK_EQUAL(json_escape_string("Formfeed\fFormfeed"), "Formfeed\\fFormfeed");
    KRATOS_CHECK_EQUAL(json_escape_string("Newline\nNewline"), "Newline\\nNewline");
    KRATOS_CHECK_EQUAL(json_escape_string("CarriageReturn\rCarriageReturn"), "CarriageReturn\\rCarriageReturn");
    KRATOS_CHECK_EQUAL(json_escape_string("Tab\tTab"), "Tab\\tTab");

    // Assert that control characters are escaped properly
    KRATOS_CHECK_EQUAL(json_escape_string("Null\0Null"), "Null\u0000Null");
}

KRATOS_TEST_CASE_IN_SUITE(CSVRowto_jsonTest, KratosExternalLibrariesFastSuite)
{
    csv::CSVRow row = make_csv_row(
        { "Col 1", "Col 2" },   // Fields
        { "A", "B" }            // Column names
    );

    KRATOS_CHECK_EQUAL(row.to_json(), "{\"A\":\"Col 1\",\"B\":\"Col 2\"}");
}

KRATOS_TEST_CASE_IN_SUITE(CSVRowto_jsonwithNumbers, KratosExternalLibrariesFastSuite)
{
    csv::CSVRow row = make_csv_row(
        { "1234.3", "234" },    // Fields
        { "A", "B"}             // Column names
    );

    KRATOS_CHECK_EQUAL(row.to_json(), "{\"A\":1234.3,\"B\":234}");
}

KRATOS_TEST_CASE_IN_SUITE(CSVRowto_jsonMixed, KratosExternalLibrariesFastSuite)
{
    csv::CSVRow row = make_csv_row(
        { "1234.3", "234", "ABCD", "AB1", "1337" },     // Fields
        { "A", "B", "C", "D", "E" }                     // Column names
    );

    // Full Row
    {
        KRATOS_CHECK_EQUAL(row.to_json(), "{\"A\":1234.3,\"B\":234,\"C\":\"ABCD\",\"D\":\"AB1\",\"E\":1337}");
    }

    // Subset
    {
        KRATOS_CHECK_EQUAL(row.to_json({ "B", "C" }), "{\"B\":234,\"C\":\"ABCD\"}");
        KRATOS_CHECK_EQUAL(row.to_json({ "B", "A" }), "{\"B\":234,\"A\":1234.3}");
    }
}

KRATOS_TEST_CASE_IN_SUITE(CSVRowto_json_arrayMixed, KratosExternalLibrariesFastSuite)
{
    csv::CSVRow row = make_csv_row(
        { "1234.3", "234", "ABCD", "AB1", "1337" },     // Fields
        { "A", "B", "C", "D", "E" }                     // Column names
    );

    // Full Row
    {
        KRATOS_CHECK_EQUAL(row.to_json_array(), "[1234.3,234,\"ABCD\",\"AB1\",1337]");
    }

    // Subset
    {
        KRATOS_CHECK_EQUAL(row.to_json_array({ "B", "C" }), "[234,\"ABCD\"]");
        KRATOS_CHECK_EQUAL(row.to_json_array({ "B", "A" }), "[234,1234.3]");
    }
}

// Reported in: https://github.com/vincentlaucsb/csv-parser/issues/68
KRATOS_TEST_CASE_IN_SUITE(CSVRowto_jsonwithWrongColumns, KratosExternalLibrariesFastSuite)
{
    std::string csv_string(R"(A,B,C,123,345,678,)");

    auto format = csv::CSVFormat();
    format.column_names({ "A", "B" });

    csv::CSVReader reader(format);
    reader.feed(csv_string);
    reader.end_feed();

    csv::CSVRow first_row;
    reader.read_row(first_row);

    // Since the column names provided were wrong, there won't be any data.
    // to_json() method should then produce an empty object instead of segfaulting.
    KRATOS_CHECK_EQUAL(first_row.to_json(), "{}");
    KRATOS_CHECK_EQUAL(first_row.to_json_array(), "[]");
}

} // namespace Testing.
} // namespace Kratos.
