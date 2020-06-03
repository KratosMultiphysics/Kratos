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

// KRATOS_TEST_CASE_IN_SUITE(guess_delim()Test-Pipe, KratosExternalLibrariesFastSuite) {
//     CSVGuessResult format = guess_format(
//         "./data/real_data/2009PowerStatus.txt");
//     KRATOS_CHECK(format.delim == '|');
//     KRATOS_CHECK(format.header_row == 0);
// }
//
// KRATOS_TEST_CASE_IN_SUITE(guess_delim()Test-Semi-Colon, KratosExternalLibrariesFastSuite) {
//     CSVGuessResult format = guess_format(
//         "./data/real_data/YEAR07_CBSA_NAC3.txt");
//     KRATOS_CHECK(format.delim == ';');
//     KRATOS_CHECK(format.header_row == 0);
// }
//
// KRATOS_TEST_CASE_IN_SUITE(guess_delim()Test-CSVwithComments, KratosExternalLibrariesFastSuite) {
//     CSVGuessResult format = guess_format(
//         "./data/fake_data/ints_comments.csv");
//     KRATOS_CHECK(format.delim == ',');
//     KRATOS_CHECK(format.header_row == 5);
// }

} // namespace Testing.
} // namespace Kratos.
