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
#include "csv.hpp"

namespace Kratos {

namespace Testing {

using namespace csv;

// // Construct a CSVRow and assert that its interface works as expected
// TEST_CASE("CSVRow Test", "[test_csv_row]") {
//     // Create a row of size 4
//     auto col_names = std::make_shared<internals::ColNames>(
//         std::vector<std::string>({ "A", "B", "C", "D" })
//         );
//
//     std::string str = "Col1"
//         "Col2"
//         "Col3"
//         "Col4";
//
//     std::vector<internals::StrBufferPos> splits = { 4, 8, 12 };
//
//     const CSVRow row(str, splits, col_names);
//
//     bool error_caught = false;
//
//     SECTION("size() Check") {
//         REQUIRE(row.size() == 4);
//     }
//
//     SECTION("operator[]") {
//         REQUIRE(row[1] == "Col2");
//         REQUIRE(row["B"] == "Col2");
//
//         REQUIRE(row[2] == "Col3");
//         REQUIRE(row["C"] == "Col3");
//     }
//
//     SECTION("operator[] Out of Bounds") {
//         try {
//             row[4].get<>();
//         }
//         catch (std::runtime_error&) {
//             error_caught = true;
//         }
//
//         REQUIRE(error_caught);
//     }
//
//     SECTION("operator[] Access Non-Existent Column") {
//         try {
//             row["Col5"].get<>();
//         }
//         catch (std::runtime_error&) {
//             error_caught = true;
//         }
//
//         REQUIRE(error_caught);
//     }
//
//     SECTION("Content Check") {
//         REQUIRE(std::vector<std::string>(row) ==
//             std::vector<std::string>({ "Col1", "Col2", "Col3", "Col4" }));
//     }
//
//     /** Allow get_sv() to be used with a const CSVField
//      *
//      *  See: https://github.com/vincentlaucsb/csv-parser/issues/86
//      *
//      */
//     SECTION("get_sv() Check") {
//         std::vector<std::string> content;
//
//         for (const auto& field : row) {
//             content.push_back(std::string(field.get_sv()));
//         }
//
//         REQUIRE(std::vector<std::string>(row) ==
//             std::vector<std::string>({ "Col1", "Col2", "Col3", "Col4" }));
//     }
// }
//
// // Integration test for CSVRow/CSVField
// TEST_CASE("CSVField operator==", "[test_csv_field_equal]") {
//     auto col_names = std::make_shared<internals::ColNames>(
//         std::vector<std::string>({ "A", "B", "C", "D" })
//         );
//
//     std::string str;
//     str += "1"
//         "2"
//         "3"
//         "3.14";
//
//     std::vector<internals::StrBufferPos> splits = { 1, 2, 3 };
//     CSVRow row(str, splits, col_names);
//
//     REQUIRE(row["A"] == 1);
//     REQUIRE(row["B"] == 2);
//     REQUIRE(row["C"] == 3);
//     REQUIRE(internals::is_equal(row["D"].get<long double>(), 3.14L));
// }

} // namespace Testing.
} // namespace Kratos.
