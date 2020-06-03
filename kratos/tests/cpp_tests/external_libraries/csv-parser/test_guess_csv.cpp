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

// //
// // guess_delim()
// //
// TEST_CASE("guess_delim() Test - Pipe", "[test_guess_pipe]") {
//     CSVGuessResult format = guess_format(
//         "./data/real_data/2009PowerStatus.txt");
//     REQUIRE(format.delim == '|');
//     REQUIRE(format.header_row == 0);
// }
//
// TEST_CASE("guess_delim() Test - Semi-Colon", "[test_guess_scolon]") {
//     CSVGuessResult format = guess_format(
//         "./data/real_data/YEAR07_CBSA_NAC3.txt");
//     REQUIRE(format.delim == ';');
//     REQUIRE(format.header_row == 0);
// }
//
// TEST_CASE("guess_delim() Test - CSV with Comments", "[test_guess_comment]") {
//     CSVGuessResult format = guess_format(
//         "./data/fake_data/ints_comments.csv");
//     REQUIRE(format.delim == ',');
//     REQUIRE(format.header_row == 5);
// }

} // namespace Testing.
} // namespace Kratos.
