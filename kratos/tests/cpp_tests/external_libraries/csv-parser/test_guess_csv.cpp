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
#include "utilities/string_utilities.h"
#include "csv-parser/include/csv.hpp"

namespace Kratos {

namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(guess_delimCSVwithComments, KratosExternalLibrariesFastSuite)
{
    const std::string working_dir = StringUtilities::ErasePartialString(__FILE__, "test_guess_csv.cpp");
    csv::CSVGuessResult format = csv::guess_format(Kratos::FilesystemExtensions::JoinPaths({working_dir,"data/fake_data/ints_comments.csv"}));
    KRATOS_CHECK_EQUAL(format.delim, ',');
    KRATOS_CHECK_EQUAL(format.header_row, 5);
}

} // namespace Testing.
} // namespace Kratos.
