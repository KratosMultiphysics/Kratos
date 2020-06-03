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
#include <stdio.h> // For remove()
#include <sstream>
#include <queue>
#include <list>

// External includes

// Project includes
#include "testing/testing.h"
#include "csv-parser/include/csv.hpp"

namespace Kratos {

namespace Testing {

using namespace csv;
using std::queue;
using std::vector;
using std::string;

// KRATOS_TEST_CASE_IN_SUITE(CSVCommaEscape, KratosExternalLibrariesFastSuite) {
//     std::string input = "Furthermore, this should be quoted.";
//     std::string correct = "\"Furthermore, this should be quoted.\"";
//
//     KRATOS_CHECK(csv_escape<>(input) == correct);
// }
//
// KRATOS_TEST_CASE_IN_SUITE(CSVQuoteEscape, KratosExternalLibrariesFastSuite) {
//     std::string input = "\"What does it mean to be RFC 4180 compliant?\" she asked.";
//     std::string correct = "\"\"\"What does it mean to be RFC 4180 compliant?\"\" she asked.\"";
//
//     KRATOS_CHECK(csv_escape<>(input) == correct);
// }
//
// KRATOS_TEST_CASE_IN_SUITE(CSVQuoteMinimal, KratosExternalLibrariesFastSuite) {
//     std::string input = "This should not be quoted";
//     KRATOS_CHECK(csv_escape<>(input) == input);
// }
//
// KRATOS_TEST_CASE_IN_SUITE(CSVQuoteAll, KratosExternalLibrariesFastSuite) {
//     std::string input = "This should be quoted";
//     std::string correct = "\"This should be quoted\"";
//     KRATOS_CHECK(csv_escape<>(input, false) == correct);
// }
//
// KRATOS_TEST_CASE_IN_SUITE(CSVtoStringstream, KratosExternalLibrariesFastSuite) {
//     std::stringstream out, correct;
//
//     // Build correct string
//     correct << "A,B,C" << std::endl << "\"1,1\",2,3" << std::endl;
//
//     queue<vector<string>> q;
//     q.push({ "A", "B", "C" });
//     q.push({ "1,1", "2", "3" });
//
//     auto writer = make_csv_writer(out);
//     for (; !q.empty(); q.pop())
//         writer.write_row(q.front());
//
//     KRATOS_CHECK(out.str() == correct.str());
// }
//
// TEMPLATE_TEST_CASE(CSV/TSVWriter-operator<<, KratosExternalLibrariesFastSuite) {
//     std::stringstream output, correct_comma, correct_tab;
//
//     // Build correct strings
//     correct_comma << "A,B,C" << std::endl << "\"1,1\",2,3" << std::endl;
//     correct_tab << "A\tB\tC" << std::endl << "1,1\t2\t3" << std::endl;
//
//     // Test input
//     auto test_row_1 = TestType({ "A", "B", "C" }),
//         test_row_2 = TestType({ "1,1", "2", "3" });
//
//     // CSV Writer
//     {
//         auto csv_writer = make_csv_writer(output);
//         csv_writer << test_row_1 << test_row_2;
//
//         KRATOS_CHECK(output.str() == correct_comma.str());
//     }
//
//     // TSV Writer
//     {
//         auto tsv_writer = make_tsv_writer(output);
//         tsv_writer << test_row_1 << test_row_2;
//
//         KRATOS_CHECK(output.str() == correct_tab.str());
//     }
// }

} // namespace Testing.
} // namespace Kratos.
