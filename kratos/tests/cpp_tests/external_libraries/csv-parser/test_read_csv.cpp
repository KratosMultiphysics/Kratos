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
#include "csv-parser/include/csv.hpp"

namespace Kratos {

namespace Testing {

using namespace csv;
using std::vector;
using std::string;

// TEST_CASE( "Test Parse Flags", "[test_parse_flags]" ) {
//     REQUIRE(internals::make_parse_flags(',', '"')[162] == internals::ParseFlags::QUOTE);
// }
//
// // Test Main Functions
// TEST_CASE( "Test Reading CSV From Direct Input", "[read_csv_direct]" ) {
//     auto rows = "A,B,C\r\n" // Header row
//                 "123,234,345\r\n"
//                 "1,2,3\r\n"
//                 "1,2,3"_csv;
//
//     // Expected Results
//     CSVRow row;
//     rows.read_row(row);
//     vector<string> first_row = {"123", "234", "345"};
//     REQUIRE( vector<string>(row) == first_row );
// }
//
// TEST_CASE("Assert UTF-8 Handling Works", "[read_utf8_direct]") {
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
//     REQUIRE(vector<string>(row) == first_row);
// }
//
// //! [Escaped Comma]
// TEST_CASE( "Test Escaped Comma", "[read_csv_comma]" ) {
//     auto rows = "A,B,C\r\n" // Header row
//                 "123,\"234,345\",456\r\n"
//                 "1,2,3\r\n"
//                 "1,2,3"_csv;
//
//     CSVRow row;
//     rows.read_row(row);
//     REQUIRE( vector<string>(row) ==
//         vector<string>({"123", "234,345", "456"}));
// }
// //! [Escaped Comma]
//
// TEST_CASE( "Test Escaped Newline", "[read_csv_newline]" ) {
//     auto rows = "A,B,C\r\n" // Header row
//                 "123,\"234\n,345\",456\r\n"
//                 "1,2,3\r\n"
//                 "1,2,3"_csv;
//
//     CSVRow row;
//     rows.read_row(row);
//     REQUIRE( vector<string>(row) ==
//         vector<string>({ "123", "234\n,345", "456" }) );
// }
//
// TEST_CASE( "Test Empty Field", "[read_empty_field]" ) {
//     // Per RFC 1480, escaped quotes should be doubled up
//     auto rows = "A,B,C\r\n" // Header row
//                 "123,\"\",456\r\n"_csv;
//
//     CSVRow row;
//     rows.read_row(row);
//     REQUIRE( vector<string>(row) ==
//         vector<string>({ "123", "", "456" }) );
// }
//
// //! [Parse Example]
// TEST_CASE( "Test Escaped Quote", "[read_csv_quote]" ) {
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
//         REQUIRE(vector<string>(row) == correct_row);
//     }
// }
// //! [Parse Example]
//
// TEST_CASE("Fragment Test", "[read_csv_fragments]") {
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
//         REQUIRE(vector<string>(row) == correct_row);
//     }
// }
//
// TEST_CASE("Test Whitespace Trimming", "[read_csv_trim]") {
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
//         REQUIRE(vector<string>(row) ==
//             vector<string>({ "123", "234\n,345", "456" }));
//         REQUIRE(row["A"] == "123");
//         REQUIRE(row["B"] == "234\n,345");
//         REQUIRE(row["C"] == "456");
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
// TEST_CASE("Test Whitespace Trimming w/ Empty Fields") {
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
//         REQUIRE(row[0].get<uint32_t>() == 1);
//         REQUIRE(row[1].get<std::string>() == "two");
//         REQUIRE(row[2].get<uint32_t>() == 3);
//
//         // Second Row
//         rows.read_row(row);
//         REQUIRE(row[0].get<uint32_t>() == 4);
//         REQUIRE(row[1].is_null());
//         REQUIRE(row[2].get<uint32_t>() == 5);
//
//         // Third Row
//         rows.read_row(row);
//         REQUIRE(row[0].is_null());
//         REQUIRE(row[1].get<uint32_t>() == 6);
//         REQUIRE(row[2].is_null());
//
//         // Fourth Row
//         rows.read_row(row);
//         REQUIRE(row[0].get<uint32_t>() == 7);
//         REQUIRE(row[1].get<uint32_t>() == 8);
//         REQUIRE(row[2].get<uint32_t>() == 9);
//     }
// }
//
// TEST_CASE("Test Variable Row Length Handling", "[read_csv_var_len]") {
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
//         REQUIRE(error_caught);
//         REQUIRE(i == 2);
//         REQUIRE(error_message.substr(0, 14) == "Line too short");
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
//         REQUIRE(rows.size() == 3);
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
//         REQUIRE(rows.size() == 5);
//         REQUIRE(rows[2][0] == 6);
//         REQUIRE(rows[2][1] == 9);
//
//         // Should be able to index extra columns via numeric index
//         REQUIRE(rows[3][2] == 7);
//         REQUIRE(rows[3][3] == 10);
//     }
// }
//
// TEST_CASE("Test read_row() CSVField - Memory", "[read_row_csvf2]") {
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
//     REQUIRE((row[0].is_float() && row[0].is_num()));
//     REQUIRE(row[0].get<std::string>().substr(0, 4) == "3.14");
//     REQUIRE(internals::is_equal(row[0].get<double>(), 3.14));
//
//     // Second Row
//     rows.read_row(row);
//     REQUIRE((row[0].is_int() && row[0].is_num()));
//     REQUIRE((row[1].is_int() && row[1].is_num()));
//     REQUIRE(row[0].get<std::string>() == "60");
//     REQUIRE(row[1].get<std::string>() == "70");
//
//     // Third Row
//     rows.read_row(row);
//     REQUIRE(row[0].is_null());
//     REQUIRE(row[1].is_null());
// }
//
// // Reported in: https://github.com/vincentlaucsb/csv-parser/issues/56
// TEST_CASE("Leading Empty Field Regression", "[empty_field_regression]") {
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
//     REQUIRE(reader.read_row(first_row));
//     REQUIRE(reader.read_row(second_row));
//
//     REQUIRE(first_row["category"] == "");
//     REQUIRE(first_row["subcategory"] == "");
//     REQUIRE(first_row["project name"] == "foo-project");
//
//     REQUIRE(second_row["category"] == "bar-category");
//     REQUIRE(second_row["subcategory"] == "");
//     REQUIRE(second_row["project name"] == "bar-project");
// }
//
// TEST_CASE("Test Parsing CSV with Dummy Column", "[read_csv_dummy]") {
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
//     REQUIRE(reader.get_col_names() == std::vector<std::string>({"A","B","C",""}));
//
//     reader.read_row(first_row);
//     REQUIRE(std::vector<std::string>(first_row) == std::vector<std::string>({
//         "123", "345", "678", ""
//     }));
// }
//
// // Reported in: https://github.com/vincentlaucsb/csv-parser/issues/67
// TEST_CASE("Comments in Header Regression", "[comments_in_header_regression]") {
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
//     REQUIRE(expected == reader.get_col_names());
// }
//
// // Reported in: https://github.com/vincentlaucsb/csv-parser/issues/92
// TEST_CASE("Long Row Test", "[long_row_regression]") {
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
//     REQUIRE(rows.get_col_names().size() == n_cols);
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
//         REQUIRE(field.get<>() == temp.str());
//         i++;
//     }
// }

} // namespace Testing.
} // namespace Kratos.
