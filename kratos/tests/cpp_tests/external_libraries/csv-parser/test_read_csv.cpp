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
#include "csv-parser/include/csv.hpp"

namespace Kratos {

namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(TestParseFlags, KratosExternalLibrariesFastSuite )
{
    KRATOS_CHECK_EQUAL(csv::internals::make_parse_flags(',', '"')[162], csv::internals::ParseFlags::QUOTE);
}

// Test Main Functions
KRATOS_TEST_CASE_IN_SUITE(TestReadingCSVFromDirectInput, KratosExternalLibrariesFastSuite )
{
    std::string csv_string = "A,B,C\r\n" // Header row
                "123,234,345\r\n"
                "1,2,3\r\n"
                "1,2,3";

    // Expected Results
    csv::CSVRow row;
    auto rows = csv::parse(csv_string);
    rows.read_row(row);
    
    std::vector<std::string> first_row = {"123", "234", "345"};
    KRATOS_CHECK_EQUAL( std::vector<std::string>(row), first_row );
}

KRATOS_TEST_CASE_IN_SUITE(AssertUTF8HandlingWorks, KratosExternalLibrariesFastSuite)
{
    // TODO: Actually check to see if flag is set
    std::string csv_string = "\uFEFFA,B,C\r\n" // Header row
        "123,234,345\r\n"
        "1,2,3\r\n"
        "1,2,3";

    // Expected Results
    csv::CSVRow row;
    auto rows = csv::parse(csv_string);
    rows.read_row(row);

    std::vector<std::string> first_row = { "123", "234", "345" };
    KRATOS_CHECK_EQUAL(std::vector<std::string>(row), first_row);
}

KRATOS_TEST_CASE_IN_SUITE(TestEscapedComma, KratosExternalLibrariesFastSuite )
{
    std::string csv_string = "A,B,C\r\n" // Header row
                "123,\"234,345\",456\r\n"
                "1,2,3\r\n"
                "1,2,3";

    csv::CSVRow row;
    auto rows = csv::parse(csv_string);
    rows.read_row(row);

    KRATOS_CHECK_EQUAL( std::vector<std::string>(row),
        std::vector<std::string>({"123", "234,345", "456"}));
}

KRATOS_TEST_CASE_IN_SUITE(TestEscapedNewline, KratosExternalLibrariesFastSuite )
{
    std::string csv_string = "A,B,C\r\n" // Header row
                "123,\"234\n,345\",456\r\n"
                "1,2,3\r\n"
                "1,2,3";

    csv::CSVRow row;
    auto rows = csv::parse(csv_string);
    rows.read_row(row);

    KRATOS_CHECK_EQUAL( std::vector<std::string>(row),
        std::vector<std::string>({ "123", "234\n,345", "456" }) );
}

KRATOS_TEST_CASE_IN_SUITE(TestEmptyField, KratosExternalLibrariesFastSuite )
{
    // Per RFC 1480, escaped quotes should be doubled up
    std::string csv_string = "A,B,C\r\n" // Header row
                "123,\"\",456\r\n";

    csv::CSVRow row;
    auto rows = csv::parse(csv_string);
    rows.read_row(row);

    KRATOS_CHECK_EQUAL( std::vector<std::string>(row),
        std::vector<std::string>({ "123", "", "456" }) );
}

KRATOS_TEST_CASE_IN_SUITE(TestEscapedQuote, KratosExternalLibrariesFastSuite )
{
    // Per RFC 1480, escaped quotes should be doubled up
    std::string csv_string = (
        "A,B,C\r\n" // Header row
        "123,\"234\"\"345\",456\r\n"
        "123,\"234\"345\",456\r\n" // Unescaped single quote (not strictly valid)
        "123,\"234\"345\",\"456\"" // Quoted field at the end
    );

    csv::CSVRow row;
    auto rows = csv::parse(csv_string);
    rows.read_row(row);

    // Expected Results: Double " is an escape for a single "
    std::vector<std::string> correct_row = {"123", "234\"345", "456"};
    for (auto& row : rows) {
        KRATOS_CHECK_EQUAL(std::vector<std::string>(row), correct_row);
    }
}

KRATOS_TEST_CASE_IN_SUITE(FragmentTest, KratosExternalLibrariesFastSuite)
{
    csv::CSVReader reader;

    reader.feed("A,B,C\r\n" // Header row
        "123,\"234\"\"345\",456\r\n");
    reader.feed("123,\"234\"345\",456\r\n"
                "123,\"234\"345\",\"456\"");
    reader.end_feed();

    // Expected Results: Double " is an escape for a single "
    std::vector<std::string> correct_row = { "123", "234\"345", "456" };
    for (auto& row : reader) {
        KRATOS_CHECK_EQUAL(std::vector<std::string>(row), correct_row);
    }
}

KRATOS_TEST_CASE_IN_SUITE(VariableRowLengthHandling, KratosExternalLibrariesFastSuite)
{
    std::string csv_string("A,B,C\r\n" // Header row
        "123,234,345\r\n"
        "1,2,3\r\n"
        "6,9\r\n" // Short row
        "6,9,7,10\r\n" // Long row
        "1,2,3"),
        error_message = "";
    bool error_caught = false;

    // Throw Error
    {
        csv::CSVFormat format;
        format.variable_columns(csv::VariableColumnPolicy::THROW);

        auto rows = parse(csv_string, format);
        size_t i = 0;

        try {
            for (auto it = rows.begin(); it != rows.end(); ++it) {
                i++;
            }
        }
        catch (std::runtime_error& err) {
            error_caught = true;
            error_message = err.what();
        }

        KRATOS_CHECK(error_caught);
        KRATOS_CHECK_EQUAL(i, 2);
        KRATOS_CHECK_EQUAL(error_message.substr(0, 14), "Line too short");
    }

    error_caught = false;

    // Ignore Row
    {
        csv::CSVFormat format;
        format.variable_columns(false);

        auto reader = parse(csv_string, format);
        std::vector<csv::CSVRow> rows(reader.begin(), reader.end());

        // Expect short/long rows to be dropped
        KRATOS_CHECK_EQUAL(rows.size(), 3);
    }

    error_caught = false;

    // Keep Row
    {
        csv::CSVFormat format;
        format.variable_columns(true);

        auto reader = parse(csv_string, format);
        std::vector<csv::CSVRow> rows(reader.begin(), reader.end());

        // Expect short/long rows to be kept
        KRATOS_CHECK_EQUAL(rows.size(), 5);
        KRATOS_CHECK_EQUAL(rows[2][0], 6);
        KRATOS_CHECK_EQUAL(rows[2][1], 9);

        // Should be able to index extra columns via numeric index
        KRATOS_CHECK_EQUAL(rows[3][2], 7);
        KRATOS_CHECK_EQUAL(rows[3][3], 10);
    }
}

KRATOS_TEST_CASE_IN_SUITE(Testread_rowCSVFieldMemory, KratosExternalLibrariesFastSuite)
{
    csv::CSVFormat format;
    format.column_names({ "A", "B" });

    std::stringstream csv_string;
    csv_string << "3.14,9999" << std::endl
        << "60,70" << std::endl
        << "," << std::endl;

    auto rows = parse(csv_string.str(), format);
    csv::CSVRow row;
    rows.read_row(row);

    // First Row
    KRATOS_CHECK((row[0].is_float() && row[0].is_num()));
    KRATOS_CHECK_EQUAL(row[0].get<std::string>().substr(0, 4), "3.14");
    KRATOS_CHECK(csv::internals::is_equal(row[0].get<double>(), 3.14));

    // Second Row
    rows.read_row(row);
    KRATOS_CHECK((row[0].is_int() && row[0].is_num()));
    KRATOS_CHECK((row[1].is_int() && row[1].is_num()));
    KRATOS_CHECK_EQUAL(row[0].get<std::string>(), "60");
    KRATOS_CHECK_EQUAL(row[1].get<std::string>(), "70");

    // Third Row
    rows.read_row(row);
    KRATOS_CHECK(row[0].is_null());
    KRATOS_CHECK(row[1].is_null());
}

// Reported in: https://github.com/vincentlaucsb/csv-parser/issues/56
KRATOS_TEST_CASE_IN_SUITE(LeadingEmptyFieldRegression, KratosExternalLibrariesFastSuite)
{
    std::string csv_string(R"(category,subcategory,project name
,,foo-project
bar-category,,bar-project
	)");
    auto format = csv::CSVFormat();
    csv::CSVReader reader(format);
    reader.feed(csv_string);
    reader.end_feed();

    csv::CSVRow first_row, second_row;
    KRATOS_CHECK(reader.read_row(first_row));
    KRATOS_CHECK(reader.read_row(second_row));

    KRATOS_CHECK_EQUAL(first_row["category"], "");
    KRATOS_CHECK_EQUAL(first_row["subcategory"], "");
    KRATOS_CHECK_EQUAL(first_row["project name"], "foo-project");

    KRATOS_CHECK_EQUAL(second_row["category"], "bar-category");
    KRATOS_CHECK_EQUAL(second_row["subcategory"], "");
    KRATOS_CHECK_EQUAL(second_row["project name"], "bar-project");
}

// Reported in: https://github.com/vincentlaucsb/csv-parser/issues/67
KRATOS_TEST_CASE_IN_SUITE(CommentsinHeaderRegression, KratosExternalLibrariesFastSuite)
{
    std::string csv_string(R"(# some extra metadata
    # some extra metadata
    timestamp,distance,angle,amplitude
    22857782,30000,-3141.59,0
    22857786,30000,-3141.09,0
    )");

    auto format = csv::CSVFormat();
    format.header_row(2);

    csv::CSVReader reader(format);
    reader.feed(csv_string);
    reader.end_feed();

    std::vector<std::string> expected = {
        "timestamp", "distance", "angle", "amplitude"
    };

    // Original issue: Leading comments appeared in column names
    const std::vector<std::string>& r_solution = reader.get_col_names();
    for (std::size_t i = 0; i < 4; ++i) {
        KRATOS_CHECK_STRING_EQUAL(expected[i], StringUtilities::RemoveWhiteSpaces(r_solution[i]));
    }
}

// Reported in: https://github.com/vincentlaucsb/csv-parser/issues/92
KRATOS_TEST_CASE_IN_SUITE(LongRowTest, KratosExternalLibrariesFastSuite)
{
    std::stringstream csv_string;
    constexpr int n_cols = 100000;

    // Make header row
    for (int i = 0; i < n_cols; i++) {
        csv_string << i;
        if (i + 1 == n_cols) {
            csv_string << std::endl;
        }
        else {
            csv_string << ',';
        }
    }

    // Make data row
    for (int i = 0; i < n_cols; i++) {
        csv_string << (double)i * 0.000001;
        if (i + 1 == n_cols) {
            csv_string << std::endl;
        }
        else {
            csv_string << ',';
        }
    }

    auto rows = csv::parse(csv_string.str());
    KRATOS_CHECK_EQUAL(rows.get_col_names().size(), n_cols);

    csv::CSVRow row;
    rows.read_row(row);

    int i = 0;

    // Make sure all CSV fields are correct
    for (auto& field : row) {
        std::stringstream temp;
        temp << (double)i * 0.000001;
        KRATOS_CHECK_EQUAL(field.get<>(), temp.str());
        i++;
    }
}

} // namespace Testing.
} // namespace Kratos.
