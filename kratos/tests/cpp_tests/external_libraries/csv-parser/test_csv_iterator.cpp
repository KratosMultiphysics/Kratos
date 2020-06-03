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
#include "csv_parser/include/csv.hpp"

namespace Kratos {

namespace Testing {

using namespace csv;

//////////////////////
// CSVRow Iterators //
//////////////////////

// TEST_CASE("Test CSVRow Interator", "[test_csv_row_iter]") {
//     auto rows = "A,B,C\r\n" // Header row
//         "123,234,345\r\n"
//         "1,2,3\r\n"
//         "1,2,3"_csv;
//
//
//     CSVRow row;
//     rows.read_row(row);
//
//     SECTION("Forwards and Backwards Iterators") {
//         // Forwards
//         REQUIRE(row.begin()->get<int>() == 123);
//         REQUIRE((row.end() - 1)->get<>() == "345");
//
//         size_t i = 0;
//         for (auto it = row.begin(); it != row.end(); ++it) {
//             if (i == 0) REQUIRE(it->get<>() == "123");
//             else if (i == 1) REQUIRE(it->get<>() == "234");
//             else  REQUIRE(it->get<>() == "345");
//
//             i++;
//         }
//
//         // Backwards
//         REQUIRE(row.rbegin()->get<int>() == 345);
//         REQUIRE((row.rend() - 1)->get<>() == "123");
//     }
//
//     SECTION("Iterator Arithmetic") {
//         REQUIRE(row.begin()->get<int>() == 123);
//         REQUIRE((row.end() - 1)->get<>() == "345");
//
//         auto row_start = row.begin();
//         REQUIRE(*(row_start + 1) == "234");
//         REQUIRE(*(row_start + 2) == "345");
//
//     }
//
//     SECTION("Post-Increment Iterator") {
//         auto it = row.begin();
//
//         REQUIRE(it++->get<int>() == 123);
//         REQUIRE(it->get<int>() == 234);
//
//         REQUIRE(it--->get<int>() == 234);
//         REQUIRE(it->get<int>() == 123);
//     }
//
//     SECTION("Range Based For") {
//         size_t i = 0;
//         for (auto& field : row) {
//             if (i == 0) REQUIRE(field.get<>() == "123");
//             else if (i == 1) REQUIRE(field.get<>() == "234");
//             else  REQUIRE(field.get<>() == "345");
//
//             i++;
//         }
//     }
// }
//
// /////////////////////////
// // CSVReader Iterators //
// /////////////////////////
//
// //! [CSVReader Iterator 1]
// TEST_CASE("Basic CSVReader Iterator Test", "[read_ints_iter]") {
//     // A file with 100 rows and columns A, B, ... J
//     // where every value in the ith row is the number i
//     CSVReader reader("./data/fake_data/ints.csv");
//     std::vector<std::string> col_names = {
//         "A", "B", "C", "D", "E", "F", "G", "H", "I", "J"
//     };
//     size_t i = 1;
//
//     SECTION("Basic Iterator") {
//         for (auto it = reader.begin(); it != reader.end(); ++it) {
//             REQUIRE((*it)[0].get<int>() == i);
//             i++;
//         }
//     }
//
//     SECTION("Iterator Post-Increment") {
//         auto it = reader.begin();
//         REQUIRE((it++)->operator[]("A").get<int>() == 1);
//         REQUIRE(it->operator[]("A").get<int>() == 2);
//     }
//
//     SECTION("Range-Based For Loop") {
//         for (auto& row : reader) {
//             for (auto& j : col_names) REQUIRE(row[j].get<int>() == i);
//             i++;
//         }
//     }
// }
// //! [CSVReader Iterator 1]
//
// //! [CSVReader Iterator 2]
// TEST_CASE("CSVReader Iterator + std::max_elem", "[iter_max_elem]") {
//     // The first is such that each value in the ith row is the number i
//     // There are 100 rows
//     // The second file is a database of California state employee salaries
//     CSVReader r1("./data/fake_data/ints.csv"),
//         r2("./data/real_data/2015_StateDepartment.csv");
//
//     // Find largest number
//     auto int_finder = [](CSVRow& left, CSVRow& right) {
//         return (left["A"].get<int>() < right["A"].get<int>());
//     };
//
//     auto max_int = std::max_element(r1.begin(), r2.end(), int_finder);
//
//     // Find highest salary
//     auto wage_finder = [](CSVRow& left, CSVRow& right) {
//         return (left["Total Wages"].get<double>() < right["Total Wages"].get<double>());
//     };
//
//     auto max_wage = std::max_element(r2.begin(), r2.end(), wage_finder);
//
//     REQUIRE((*max_int)["A"] == 100);
//     REQUIRE((*max_wage)["Total Wages"] == "812064.87");
// }
// //! [CSVReader Iterator 2]

} // namespace Testing.
} // namespace Kratos.
