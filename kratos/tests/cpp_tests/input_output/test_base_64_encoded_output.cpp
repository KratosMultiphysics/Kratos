//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <sstream>
#include <vector>

// External includes

// Project includes
#include "testing/testing.h"
#include "input_output/base_64_encoded_output.h"

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(Base64EcodedEmptyOutput, KratosCoreFastSuite)
{
    std::stringstream output;

    {
        std::vector<uint8_t> data = {};
        auto encoder = Base64EncodedOutput(output);
        encoder.WriteData(data.begin(), data.size());
    }

    KRATOS_CHECK_EQUAL(output.str(), "");
}

KRATOS_TEST_CASE_IN_SUITE(Base64EcodedSingleCharacter, KratosCoreFastSuite)
{
    std::stringstream output;

    {
        std::vector<uint8_t> data = {'A'};
        auto encoder = Base64EncodedOutput(output);
        encoder.WriteData(data.begin(), data.size());
    }

    KRATOS_CHECK_EQUAL(output.str(), "QQ==");
}

KRATOS_TEST_CASE_IN_SUITE(Base64EcodedExampleString, KratosCoreFastSuite)
{
    std::stringstream output;

    {
        std::string data = "Man";
        auto encoder = Base64EncodedOutput(output);
        encoder.WriteData(data.begin(), data.size());
    }

    KRATOS_CHECK_EQUAL(output.str(), "TWFu");
}

KRATOS_TEST_CASE_IN_SUITE(Base64EcodedPadding, KratosCoreFastSuite)
{
    std::stringstream output;

    {
        std::vector<uint8_t> data = {'A', 'B', 'C'};
        auto encoder = Base64EncodedOutput(output);
        encoder.WriteData(data.begin(), data.size());
    }

    KRATOS_CHECK_EQUAL(output.str(), "QUJD");
}

KRATOS_TEST_CASE_IN_SUITE(Base64EcodedLargeInput, KratosCoreFastSuite)
{
    std::stringstream output;

    {
        std::vector<uint8_t> data(1000, 'U');
        auto encoder = Base64EncodedOutput(output);
        encoder.WriteData(data.begin(), data.size());
    }

    std::string expected(1000 / 3 * 4, 'V');
    expected += "VQ==";
    KRATOS_CHECK_EQUAL(output.str(), expected);
}

KRATOS_TEST_CASE_IN_SUITE(Base64EcodedSmallInput, KratosCoreFastSuite)
{
    std::stringstream output;

    {
        auto encoder = Base64EncodedOutput(output);

        const unsigned int v1 = 1;
        encoder.WriteData(&v1, 1);

        const char v2 = 1;
        encoder.WriteData(&v2, 1);
    }

    KRATOS_CHECK_EQUAL(output.str(), "AQAAAAE=");
}

} // namespace Kratos::Testing