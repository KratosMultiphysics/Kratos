//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes
#include <numeric>
#include <vector>

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "utilities/vectorization_helper.h"

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(VectorizationHelperScalarOutput, KratosCoreFastSuite)
{
    std::vector<double> container(100);
    std::iota(container.begin(), container.end(), 0.0);

    const double offset = 2.5;

    auto p_result = VectorizationHelper(
        container, [](const double rValue, const double Offset) { return rValue * 3.0 + Offset; }, offset);

    const auto& r_shape = p_result->Shape();
    KRATOS_EXPECT_EQ(r_shape.size(), 1u);
    KRATOS_EXPECT_EQ(r_shape[0], 100u);

    const auto data = p_result->ViewData();
    for (std::size_t i = 0; i < container.size(); ++i) {
        KRATOS_EXPECT_DOUBLE_EQ(data[i], container[i] * 3.0 + offset);
    }
}

KRATOS_TEST_CASE_IN_SUITE(VectorizationHelperFixedArrayOutput, KratosCoreFastSuite)
{
    std::vector<double> container(50);
    std::iota(container.begin(), container.end(), 0.0);

    const double scale = 2.0;

    auto p_result = VectorizationHelper(
        container,
        [](const double& rValue, const double Scale) {
            return array_1d<double, 3>({rValue * Scale, rValue * Scale, rValue * Scale});
        },
        scale);

    const auto& r_shape = p_result->Shape();
    KRATOS_EXPECT_EQ(r_shape.size(), 2u);
    KRATOS_EXPECT_EQ(r_shape[0], 50u);
    KRATOS_EXPECT_EQ(r_shape[1], 3u);

    const auto data = p_result->ViewData();
    for (std::size_t i = 0; i < container.size(); ++i) {
        for (std::size_t j = 0; j < 3u; ++j) {
            KRATOS_EXPECT_DOUBLE_EQ(data[i * 3u + j], container[i] * scale);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(VectorizationHelperPerfectForwarding, KratosCoreFastSuite)
{
    // non-copyable argument, passed as rvalue: it must be moved once into the helper
    // and used (by reference) for every one of the calls
    struct NonCopyableArgument
    {
        explicit NonCopyableArgument(const double Value) : mValue(Value) {}
        NonCopyableArgument(const NonCopyableArgument&) = delete;
        NonCopyableArgument(NonCopyableArgument&&) = default;

        double mValue;
    };

    std::vector<double> container(37);
    std::iota(container.begin(), container.end(), 0.0);

    const double lvalue_offset = 1.0;

    auto p_result = VectorizationHelper(
        container,
        [](const double& rValue, const double Offset, const NonCopyableArgument& rNonCopyable) {
            return rValue + Offset + rNonCopyable.mValue;
        },
        lvalue_offset,
        NonCopyableArgument(100.0));

    const auto data = p_result->ViewData();
    for (std::size_t i = 0; i < container.size(); ++i) {
        KRATOS_EXPECT_DOUBLE_EQ(data[i], container[i] + 1.0 + 100.0);
    }
}

KRATOS_TEST_CASE_IN_SUITE(VectorizationHelperModelPartNodes, KratosCoreFastSuite)
{
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("TestModelPart");

    for (std::size_t i = 1; i < 21u; ++i) {
        r_model_part.CreateNewNode(static_cast<IndexType>(i), static_cast<double>(i), 0.0, 0.0);
    }

    const double scale = 4.0;

    auto p_result = VectorizationHelper(
        r_model_part.Nodes(),
        [](Node& rNode, const double Scale) { return rNode.X() * Scale; },
        scale);

    const auto& r_shape = p_result->Shape();
    KRATOS_EXPECT_EQ(r_shape.size(), 1u);
    KRATOS_EXPECT_EQ(r_shape[0], 20u);

    const auto data = p_result->ViewData();
    for (std::size_t i = 0; i < 20u; ++i) {
        KRATOS_EXPECT_DOUBLE_EQ(data[i], static_cast<double>(i + 1) * scale);
    }
}

} // namespace Kratos::Testing
