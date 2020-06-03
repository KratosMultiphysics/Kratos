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

static std::string err_preamble = "There should be no overlap between "
    "the quote character, the set of possible "
    "delimiters and the set of whitespace characters.";

// // Assert that an error is thrown if whitespace, delimiter, and quote
// KRATOS_TEST_CASE_IN_SUITE(CSVFormat-OverlappingCharacters, KratosExternalLibrariesFastSuite) {
//     CSVFormat format;
//     bool err_caught = false;
//
//     // Tab
//     {
//         try {
//             format.delimiter('\t').quote('"').trim({ '\t' });
//         }
//         catch (std::runtime_error& err) {
//             err_caught = true;
//             KRATOS_CHECK(err.what() == std::string(err_preamble + " Offending characters: '\t'."));
//         }
//
//         KRATOS_CHECK(err_caught);
//     }
//
//     // Tab with multiple other characters 
//     {
//         try {
//             format.delimiter({ ',', '\t' }).quote('"').trim({ ' ', '\t' });
//         }
//         catch (std::runtime_error& err) {
//             err_caught = true;
//             KRATOS_CHECK(err.what() == std::string(err_preamble + " Offending characters: '\t'."));
//         }
//
//         KRATOS_CHECK(err_caught);
//     }
//
//     // Repeated quote
//     {
//         try {
//             format.delimiter({ ',', '"' }).quote('"').trim({ ' ', '\t' });
//         }
//         catch (std::runtime_error& err) {
//             err_caught = true;
//             KRATOS_CHECK(err.what() == std::string(err_preamble + " Offending characters: '\"'."));
//         }
//
//         KRATOS_CHECK(err_caught);
//     }
//
//     // Multiple offenders
//     {
//         try {
//             format.delimiter({ ',', '\t', ' ' }).quote('"').trim({ ' ', '\t' });
//         }
//         catch (std::runtime_error& err) {
//             err_caught = true;
//             KRATOS_CHECK(err.what() == std::string(err_preamble + " Offending characters: '\t', ' '."));
//         }
//
//         KRATOS_CHECK(err_caught);
//     }
// }

} // namespace Testing.
} // namespace Kratos.
