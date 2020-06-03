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

// TEST_CASE("col_pos() Test", "[test_col_pos]") {
//     int pos = get_col_pos(
//         "./tests/data/real_data/2015_StateDepartment.csv",
//         "Entity Type");
//     REQUIRE(pos == 1);
// }
//
// TEST_CASE("Prevent Column Names From Being Overwritten", "[csv_col_names_overwrite]") {
//     std::vector<std::string> column_names = { "A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10" };
//
//     // Test against a variety of different CSVFormat objects
//     std::vector<CSVFormat> formats = {};
//     formats.push_back(CSVFormat::guess_csv());
//     formats.push_back(CSVFormat());
//     formats.back().delimiter(std::vector<char>({ ',', '\t', '|'}));
//     formats.push_back(CSVFormat());
//     formats.back().delimiter(std::vector<char>({ ',', '~'}));
//
//     for (auto& format_in : formats) {
//         // Set up the CSVReader
//         format_in.column_names(column_names);
//         CSVReader reader("./tests/data/fake_data/ints_comments.csv", format_in);
//
//         // Assert that column names weren't overwritten
//         CSVFormat format_out = reader.get_format();
//         REQUIRE(reader.get_col_names() == column_names);
//         REQUIRE(format_out.get_delim() == ',');
//         REQUIRE(format_out.get_header() == 5);
//     }
// }
//
// // get_file_info()
// TEST_CASE("get_file_info() Test", "[test_file_info]") {
//     CSVFileInfo info = get_file_info(
//         "./tests/data/real_data/2009PowerStatus.txt");
//
//     REQUIRE(info.delim == '|');
//     REQUIRE(info.n_rows == 37960); // Can confirm with Excel
//     REQUIRE(info.n_cols == 3);
//     REQUIRE(info.col_names == vector<string>({"ReportDt", "Unit", "Power"}));
// }
//
// TEST_CASE("Non-Existent CSV", "[read_ghost_csv]") {
//     // Make sure attempting to parse a non-existent CSV throws an error
//     bool error_caught = false;
//
//     try {
//         CSVReader reader("./lochness.csv");
//     }
//     catch (std::runtime_error& err) {
//         error_caught = true;
//         REQUIRE(err.what() == std::string("Cannot open file ./lochness.csv"));
//     }
//
//     REQUIRE(error_caught);
// }
//
// TEST_CASE( "Test Read CSV with Header Row", "[read_csv_header]" ) {
//     // Header on first row
//     const std::string data_file = "./tests/data/real_data/2015_StateDepartment.csv";
//     CSVReader reader(data_file, CSVFormat());
//     CSVRow row;
//     reader.read_row(row); // Populate row with first line
//
//     // Expected Results
//     vector<string> col_names = {
//         "Year", "Entity Type", "Entity Group", "Entity Name",
//         "Department / Subdivision", "Position", "Elected Official",
//         "Judicial", "Other Positions", "Min Classification Salary",
//         "Max Classification Salary", "Reported Base Wage", "Regular Pay",
//         "Overtime Pay", "Lump-Sum Pay", "Other Pay", "Total Wages",
//         "Defined Benefit Plan Contribution", "Employees Retirement Cost Covered",
//         "Deferred Compensation Plan", "Health Dental Vision",
//         "Total Retirement and Health Cost", "Pension Formula",
//         "Entity URL", "Entity Population", "Last Updated",
//         "Entity County", "Special District Activities"
//     };
//
//     vector<string> first_row = {
//         "2015","State Department","","Administrative Law, Office of","",
//         "Assistant Chief Counsel","False","False","","112044","129780",""
//         ,"133020.06","0","2551.59","2434.8","138006.45","34128.65","0","0"
//         ,"15273.97","49402.62","2.00% @ 55","http://www.spb.ca.gov/","",
//         "08/02/2016","",""
//     };
//
//     REQUIRE( vector<string>(row) == first_row );
//     REQUIRE( get_col_names(data_file) == col_names );
//
//     // Skip to end
//     while (reader.read_row(row));
//     REQUIRE( reader.num_rows == 246497 );
// }
//
// //
// // read_row()
// //
// //! [CSVField Example]
// TEST_CASE("Test read_row() CSVField - Easy", "[read_row_csvf1]") {
//     // Test that integers are type-casted properly
//     CSVReader reader("./tests/data/fake_data/ints.csv");
//     CSVRow row;
//
//     while (reader.read_row(row)) {
//         for (size_t i = 0; i < row.size(); i++) {
//             REQUIRE(row[i].is_int());
//             REQUIRE(row[i].get<int>() <= 100);
//         }
//     }
// }
// //! [CSVField Example]
//
// TEST_CASE("Test read_row() CSVField - Power Status", "[read_row_csvf3]") {
//     CSVReader reader("./tests/data/real_data/2009PowerStatus.txt");
//     CSVRow row;
//
//     size_t date = reader.index_of("ReportDt"),
//         unit = reader.index_of("Unit"),
//         power = reader.index_of("Power");
//
//     // Try to find a non-existent column
//     REQUIRE(reader.index_of("metallica") == CSV_NOT_FOUND);
//
//     for (size_t i = 0; reader.read_row(row); i++) {
//         // Assert correct types
//         REQUIRE(row[date].is_str());
//         REQUIRE(row[unit].is_str());
//         REQUIRE(row[power].is_int());
//
//         // Spot check
//         if (i == 2) {
//             REQUIRE(row[power].get<int>() == 100);
//             REQUIRE(row[date].get<>() == "12/31/2009"); // string_view
//             REQUIRE(row[unit].get<std::string>() == "Beaver Valley 1");
//         }
//     }
// }

} // namespace Testing.
} // namespace Kratos.
