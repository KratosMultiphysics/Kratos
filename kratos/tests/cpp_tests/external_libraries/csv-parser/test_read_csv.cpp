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
#include "csv_parser/include/csv.hpp"

namespace Kratos {

namespace Testing {

using namespace csv;
using std::vector;
using std::string;

// KRATOS_TEST_CASE_IN_SUITE( "Test Parse Flags", "[test_parse_flags]" ) {
//     KRATOS_CHECK(internals::make_parse_flags(',', '"')[162] == internals::ParseFlags::QUOTE);
// }
//
// // Test Main Functions
// KRATOS_TEST_CASE_IN_SUITE( "Test Reading CSV From Direct Input", "[read_csv_direct]" ) {
//     auto rows = "A,B,C\r\n" // Header row
//                 "123,234,345\r\n"
//                 "1,2,3\r\n"
//                 "1,2,3"_csv;
//
//     // Expected Results
//     CSVRow row;
//     rows.read_row(row);
//     vector<string> first_row = {"123", "234", "345"};
//     KRATOS_CHECK( vector<string>(row) == first_row );
// }
//
// KRATOS_TEST_CASE_IN_SUITE("Assert UTF-8 Handling Works", "[read_utf8_direct]") {
//     // TODO: Actually check to see if flag is set
//     auto rows = "\uFEFFA,B,C\r\n" // Header row
//         "123,234,345\r\n"
//         "1,2,3\r\n"
//         "1,2,3"_csv;
//
//     // Expected Results
//     CSVRow row;
//     rows.read_row(row);
//     vector<string> first_row = { "123", "234", "345" };
//     KRATOS_CHECK(vector<string>(row) == first_row);
// }
//
// //! [Escaped Comma]
// KRATOS_TEST_CASE_IN_SUITE( "Test Escaped Comma", "[read_csv_comma]" ) {
//     auto rows = "A,B,C\r\n" // Header row
//                 "123,\"234,345\",456\r\n"
//                 "1,2,3\r\n"
//                 "1,2,3"_csv;
//
//     CSVRow row;
//     rows.read_row(row);
//     KRATOS_CHECK( vector<string>(row) ==
//         vector<string>({"123", "234,345", "456"}));
// }
// //! [Escaped Comma]
//
// KRATOS_TEST_CASE_IN_SUITE( "Test Escaped Newline", "[read_csv_newline]" ) {
//     auto rows = "A,B,C\r\n" // Header row
//                 "123,\"234\n,345\",456\r\n"
//                 "1,2,3\r\n"
//                 "1,2,3"_csv;
//
//     CSVRow row;
//     rows.read_row(row);
//     KRATOS_CHECK( vector<string>(row) ==
//         vector<string>({ "123", "234\n,345", "456" }) );
// }
//
// KRATOS_TEST_CASE_IN_SUITE( "Test Empty Field", "[read_empty_field]" ) {
//     // Per RFC 1480, escaped quotes should be doubled up
//     auto rows = "A,B,C\r\n" // Header row
//                 "123,\"\",456\r\n"_csv;
//
//     CSVRow row;
//     rows.read_row(row);
//     KRATOS_CHECK( vector<string>(row) ==
//         vector<string>({ "123", "", "456" }) );
// }
//
// //! [Parse Example]
// KRATOS_TEST_CASE_IN_SUITE( "Test Escaped Quote", "[read_csv_quote]" ) {
//     // Per RFC 1480, escaped quotes should be doubled up
//     string csv_string = (
//         "A,B,C\r\n" // Header row
//         "123,\"234\"\"345\",456\r\n"
//         "123,\"234\"345\",456\r\n" // Unescaped single quote (not strictly valid)
//         "123,\"234\"345\",\"456\"" // Quoted field at the end
//     );
//
//     auto rows = parse(csv_string);
//
//     // Expected Results: Double " is an escape for a single "
//     vector<string> correct_row = {"123", "234\"345", "456"};
//     for (auto& row : rows) {
//         KRATOS_CHECK(vector<string>(row) == correct_row);
//     }
// }
// //! [Parse Example]
//
// KRATOS_TEST_CASE_IN_SUITE("Fragment Test", "[read_csv_fragments]") {
//     CSVReader reader;
//
//     reader.feed("A,B,C\r\n" // Header row
//         "123,\"234\"\"345\",456\r\n");
//     reader.feed("123,\"234\"345\",456\r\n"
//                 "123,\"234\"345\",\"456\"");
//     reader.end_feed();
//
//     // Expected Results: Double " is an escape for a single "
//     vector<string> correct_row = { "123", "234\"345", "456" };
//     for (auto& row : reader) {
//         KRATOS_CHECK(vector<string>(row) == correct_row);
//     }
// }
//
// KRATOS_TEST_CASE_IN_SUITE("Test Whitespace Trimming", "[read_csv_trim]") {
//     auto row_str = GENERATE(as<std::string> {},
//         "A,B,C\r\n" // Header row
//         "123,\"234\n,345\",456\r\n",
//
//         // Random spaces
//         "A,B,C\r\n"
//         "   123,\"234\n,345\",    456\r\n",
//
//         // Random spaces + tabs
//         "A,B,C\r\n"
//         "\t\t   123,\"234\n,345\",    456\r\n",
//
//         // Spaces in quote escaped field
//         "A,B,C\r\n"
//         "\t\t   123,\"   234\n,345  \t\",    456\r\n",
//
//         // Spaces in one header column
//         "A,B,        C\r\n"
//         "123,\"234\n,345\",456\r\n",
//
//         // Random spaces + tabs in header
//         "\t A,  B\t,     C\r\n"
//         "123,\"234\n,345\",456\r\n",
//
//         // Random spaces in header + data
//         "A,B,        C\r\n"
//         "123,\"234\n,345\",  456\r\n"
//     );
//
//     SECTION("Parse Test") {
//         CSVFormat format;
//         format.header_row(0)
//             .trim({ '\t', ' ' })
//             .delimiter(',');
//
//         auto rows = parse(row_str, format);
//         CSVRow row;
//         rows.read_row(row);
//
//         KRATOS_CHECK(vector<string>(row) ==
//             vector<string>({ "123", "234\n,345", "456" }));
//         KRATOS_CHECK(row["A"] == "123");
//         KRATOS_CHECK(row["B"] == "234\n,345");
//         KRATOS_CHECK(row["C"] == "456");
//     }
// }
//
// std::vector<std::string> make_whitespace_test_cases() {
//     std::vector<std::string> test_cases = {};
//     std::stringstream ss;
//
//     ss << "1, two,3" << std::endl
//         << "4, ,5" << std::endl
//         << " ,6, " << std::endl
//         << "7,8,9 " << std::endl;
//     test_cases.push_back(ss.str());
//     ss.clear();
//
//     // Lots of Whitespace
//     ss << "1, two,3" << std::endl
//         << "4,                    ,5" << std::endl
//         << "         ,6,       " << std::endl
//         << "7,8,9 " << std::endl;
//     test_cases.push_back(ss.str());
//     ss.clear();
//
//     // Same as above but there's whitespace around 6
//     ss << "1, two,3" << std::endl
//         << "4,                    ,5" << std::endl
//         << "         , 6 ,       " << std::endl
//         << "7,8,9 " << std::endl;
//     test_cases.push_back(ss.str());
//     ss.clear();
//
//     // Tabs
//     ss << "1, two,3" << std::endl
//         << "4, \t ,5" << std::endl
//         << "\t\t\t\t\t ,6, \t " << std::endl
//         << "7,8,9 " << std::endl;
//     test_cases.push_back(ss.str());
//     ss.clear();
//
//     return test_cases;
// }
//
// KRATOS_TEST_CASE_IN_SUITE("Test Whitespace Trimming w/ Empty Fields") {
//     auto csv_string = GENERATE(from_range(make_whitespace_test_cases()));
//
//     SECTION("Parse Test") {
//         CSVFormat format;
//         format.column_names({ "A", "B", "C" })
//             .trim({ ' ', '\t' });
//
//         auto rows = parse(csv_string, format);
//         CSVRow row;
//
//         // First Row
//         rows.read_row(row);
//         KRATOS_CHECK(row[0].get<uint32_t>() == 1);
//         KRATOS_CHECK(row[1].get<std::string>() == "two");
//         KRATOS_CHECK(row[2].get<uint32_t>() == 3);
//
//         // Second Row
//         rows.read_row(row);
//         KRATOS_CHECK(row[0].get<uint32_t>() == 4);
//         KRATOS_CHECK(row[1].is_null());
//         KRATOS_CHECK(row[2].get<uint32_t>() == 5);
//
//         // Third Row
//         rows.read_row(row);
//         KRATOS_CHECK(row[0].is_null());
//         KRATOS_CHECK(row[1].get<uint32_t>() == 6);
//         KRATOS_CHECK(row[2].is_null());
//
//         // Fourth Row
//         rows.read_row(row);
//         KRATOS_CHECK(row[0].get<uint32_t>() == 7);
//         KRATOS_CHECK(row[1].get<uint32_t>() == 8);
//         KRATOS_CHECK(row[2].get<uint32_t>() == 9);
//     }
// }
//
// KRATOS_TEST_CASE_IN_SUITE("Test Variable Row Length Handling", "[read_csv_var_len]") {
//     string csv_string("A,B,C\r\n" // Header row
//         "123,234,345\r\n"
//         "1,2,3\r\n"
//         "6,9\r\n" // Short row
//         "6,9,7,10\r\n" // Long row
//         "1,2,3"),
//         error_message = "";
//     bool error_caught = false;
//
//     SECTION("Throw Error") {
//         CSVFormat format;
//         format.variable_columns(VariableColumnPolicy::THROW);
//
//         auto rows = parse(csv_string, format);
//         size_t i = 0;
//
//         try {
//             for (auto it = rows.begin(); it != rows.end(); ++it) {
//                 i++;
//             }
//         }
//         catch (std::runtime_error& err) {
//             error_caught = true;
//             error_message = err.what();
//         }
//
//         KRATOS_CHECK(error_caught);
//         KRATOS_CHECK(i == 2);
//         KRATOS_CHECK(error_message.substr(0, 14) == "Line too short");
//     }
//
//     SECTION("Ignore Row") {
//         CSVFormat format;
//         format.variable_columns(false);
//
//         auto reader = parse(csv_string, format);
//         std::vector<CSVRow> rows(reader.begin(), reader.end());
//
//         // Expect short/long rows to be dropped
//         KRATOS_CHECK(rows.size() == 3);
//     }
//
//     SECTION("Keep Row") {
//         CSVFormat format;
//         format.variable_columns(true);
//
//         auto reader = parse(csv_string, format);
//         std::vector<CSVRow> rows(reader.begin(), reader.end());
//
//         // Expect short/long rows to be kept
//         KRATOS_CHECK(rows.size() == 5);
//         KRATOS_CHECK(rows[2][0] == 6);
//         KRATOS_CHECK(rows[2][1] == 9);
//
//         // Should be able to index extra columns via numeric index
//         KRATOS_CHECK(rows[3][2] == 7);
//         KRATOS_CHECK(rows[3][3] == 10);
//     }
// }
//
// KRATOS_TEST_CASE_IN_SUITE("Test read_row() CSVField - Memory", "[read_row_csvf2]") {
//     CSVFormat format;
//     format.column_names({ "A", "B" });
//
//     std::stringstream csv_string;
//     csv_string << "3.14,9999" << std::endl
//         << "60,70" << std::endl
//         << "," << std::endl;
//
//     auto rows = parse(csv_string.str(), format);
//     CSVRow row;
//     rows.read_row(row);
//
//     // First Row
//     KRATOS_CHECK((row[0].is_float() && row[0].is_num()));
//     KRATOS_CHECK(row[0].get<std::string>().substr(0, 4) == "3.14");
//     KRATOS_CHECK(internals::is_equal(row[0].get<double>(), 3.14));
//
//     // Second Row
//     rows.read_row(row);
//     KRATOS_CHECK((row[0].is_int() && row[0].is_num()));
//     KRATOS_CHECK((row[1].is_int() && row[1].is_num()));
//     KRATOS_CHECK(row[0].get<std::string>() == "60");
//     KRATOS_CHECK(row[1].get<std::string>() == "70");
//
//     // Third Row
//     rows.read_row(row);
//     KRATOS_CHECK(row[0].is_null());
//     KRATOS_CHECK(row[1].is_null());
// }
//
// // Reported in: https://github.com/vincentlaucsb/csv-parser/issues/56
// KRATOS_TEST_CASE_IN_SUITE("Leading Empty Field Regression", "[empty_field_regression]") {
//     std::string csv_string(R"(category,subcategory,project name
// ,,foo-project
// bar-category,,bar-project
// 	)");
//     auto format = csv::CSVFormat();
//     csv::CSVReader reader(format);
//     reader.feed(csv_string);
//     reader.end_feed();
//
//     CSVRow first_row, second_row;
//     KRATOS_CHECK(reader.read_row(first_row));
//     KRATOS_CHECK(reader.read_row(second_row));
//
//     KRATOS_CHECK(first_row["category"] == "");
//     KRATOS_CHECK(first_row["subcategory"] == "");
//     KRATOS_CHECK(first_row["project name"] == "foo-project");
//
//     KRATOS_CHECK(second_row["category"] == "bar-category");
//     KRATOS_CHECK(second_row["subcategory"] == "");
//     KRATOS_CHECK(second_row["project name"] == "bar-project");
// }
//
// KRATOS_TEST_CASE_IN_SUITE("Test Parsing CSV with Dummy Column", "[read_csv_dummy]") {
//     std::string csv_string(R"(A,B,C,
// 123,345,678,)");
//
//     auto format = csv::CSVFormat();
//     csv::CSVReader reader(format);
//     reader.feed(csv_string);
//     reader.end_feed();
//
//     CSVRow first_row;
//
//     KRATOS_CHECK(reader.get_col_names() == std::vector<std::string>({"A","B","C",""}));
//
//     reader.read_row(first_row);
//     KRATOS_CHECK(std::vector<std::string>(first_row) == std::vector<std::string>({
//         "123", "345", "678", ""
//     }));
// }
//
// // Reported in: https://github.com/vincentlaucsb/csv-parser/issues/67
// KRATOS_TEST_CASE_IN_SUITE("Comments in Header Regression", "[comments_in_header_regression]") {
//     std::string csv_string(R"(# some extra metadata
// # some extra metadata
// timestamp,distance,angle,amplitude
// 22857782,30000,-3141.59,0
// 22857786,30000,-3141.09,0
// )");
//
//     auto format = csv::CSVFormat();
//     format.header_row(2);
//
//     csv::CSVReader reader(format);
//     reader.feed(csv_string);
//     reader.end_feed();
//
//     std::vector<std::string> expected = {
//         "timestamp", "distance", "angle", "amplitude"
//     };
//
//     // Original issue: Leading comments appeared in column names
//     KRATOS_CHECK(expected == reader.get_col_names());
// }
//
// // Reported in: https://github.com/vincentlaucsb/csv-parser/issues/92
// KRATOS_TEST_CASE_IN_SUITE("Long Row Test", "[long_row_regression]") {
//     std::stringstream csv_string;
//     constexpr int n_cols = 100000;
//
//     // Make header row
//     for (int i = 0; i < n_cols; i++) {
//         csv_string << i;
//         if (i + 1 == n_cols) {
//             csv_string << std::endl;
//         }
//         else {
//             csv_string << ',';
//         }
//     }
//
//     // Make data row
//     for (int i = 0; i < n_cols; i++) {
//         csv_string << (double)i * 0.000001;
//         if (i + 1 == n_cols) {
//             csv_string << std::endl;
//         }
//         else {
//             csv_string << ',';
//         }
//     }
//
//     auto rows = parse(csv_string.str());
//     KRATOS_CHECK(rows.get_col_names().size() == n_cols);
//
//     CSVRow row;
//     rows.read_row(row);
//
//     int i = 0;
//
//     // Make sure all CSV fields are correct
//     for (auto& field : row) {
//         std::stringstream temp;
//         temp << (double)i * 0.000001;
//         KRATOS_CHECK(field.get<>() == temp.str());
//         i++;
//     }
// }

} // namespace Testing.
} // namespace Kratos.
