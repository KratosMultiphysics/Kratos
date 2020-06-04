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
#include "includes/kratos_filesystem.h"
#include "csv-parser/include/csv.hpp"

namespace Kratos {

namespace Testing {

using namespace csv;
using std::vector;
using std::string;

/**
 * Erase First Occurrence of given  substring from main string.
 */
std::string ErasePartialString(const std::string& rMainStrint, const std::string& rToErase)
{
    // Value to return
    std::string sub_string = rMainStrint;

    // Search for the substring in string
    size_t pos = sub_string.find(rToErase);

    if (pos != std::string::npos) {
        // If found then erase it from string
        sub_string.erase(pos, rToErase.length());
    }

    return sub_string;
}

KRATOS_TEST_CASE_IN_SUITE(guess_delimPipe, KratosExternalLibrariesFastSuite)
{
    const std::string working_dir = ErasePartialString(__FILE__, "test_guess_csv.cpp");
    CSVGuessResult format = guess_format(Kratos::FilesystemExtensions::JoinPaths({working_dir,"data/real_data/2009PowerStatus.txt"}));
    KRATOS_CHECK_EQUAL(format.delim, '|');
    KRATOS_CHECK_EQUAL(format.header_row, 0);
}

KRATOS_TEST_CASE_IN_SUITE(guess_delimSemiColon, KratosExternalLibrariesFastSuite)
{
    const std::string working_dir = ErasePartialString(__FILE__, "test_guess_csv.cpp");
    CSVGuessResult format = guess_format(Kratos::FilesystemExtensions::JoinPaths({working_dir,"data/real_data/YEAR07_CBSA_NAC3.txt"}));
    KRATOS_CHECK_EQUAL(format.delim, ';');
    KRATOS_CHECK_EQUAL(format.header_row, 0);
}

KRATOS_TEST_CASE_IN_SUITE(guess_delimCSVwithComments, KratosExternalLibrariesFastSuite)
{
    const std::string working_dir = ErasePartialString(__FILE__, "test_guess_csv.cpp");
    CSVGuessResult format = guess_format(Kratos::FilesystemExtensions::JoinPaths({working_dir,"data/fake_data/ints_comments.csv"}));
    KRATOS_CHECK_EQUAL(format.delim, ',');
    KRATOS_CHECK_EQUAL(format.header_row, 5);
}

} // namespace Testing.
} // namespace Kratos.
