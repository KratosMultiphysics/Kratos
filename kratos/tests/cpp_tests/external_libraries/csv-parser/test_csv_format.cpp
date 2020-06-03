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

using namespace csv;

static std::string err_preamble = "There should be no overlap between "
    "the quote character, the set of possible "
    "delimiters and the set of whitespace characters.";

// // Assert that an error is thrown if whitespace, delimiter, and quote
// TEST_CASE("CSVFormat - Overlapping Characters", "[csv_format_overlap]") {
//     CSVFormat format;
//     bool err_caught = false;
//
//     SECTION("Tab") {
//         try {
//             format.delimiter('\t').quote('"').trim({ '\t' });
//         }
//         catch (std::runtime_error& err) {
//             err_caught = true;
//             REQUIRE(err.what() == std::string(err_preamble + " Offending characters: '\t'."));
//         }
//
//         REQUIRE(err_caught);
//     }
//
//     SECTION("Tab with multiple other characters") {
//         try {
//             format.delimiter({ ',', '\t' }).quote('"').trim({ ' ', '\t' });
//         }
//         catch (std::runtime_error& err) {
//             err_caught = true;
//             REQUIRE(err.what() == std::string(err_preamble + " Offending characters: '\t'."));
//         }
//
//         REQUIRE(err_caught);
//     }
//
//     SECTION("Repeated quote") {
//         try {
//             format.delimiter({ ',', '"' }).quote('"').trim({ ' ', '\t' });
//         }
//         catch (std::runtime_error& err) {
//             err_caught = true;
//             REQUIRE(err.what() == std::string(err_preamble + " Offending characters: '\"'."));
//         }
//
//         REQUIRE(err_caught);
//     }
//
//     SECTION("Multiple offenders") {
//         try {
//             format.delimiter({ ',', '\t', ' ' }).quote('"').trim({ ' ', '\t' });
//         }
//         catch (std::runtime_error& err) {
//             err_caught = true;
//             REQUIRE(err.what() == std::string(err_preamble + " Offending characters: '\t', ' '."));
//         }
//
//         REQUIRE(err_caught);
//     }
// }

} // namespace Testing.
} // namespace Kratos.
