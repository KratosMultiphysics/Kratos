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
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "expression/literal_flat_expression.h"
#include "expression/binary_expression.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(ExpressionIteratorLazy, KratosCoreFastSuite) {
    auto p_literal_scalar_expression = LiteralFlatExpression<double>::Create(10, {});
    auto p_literal_array_3_expression = LiteralFlatExpression<double>::Create(10, {3});

    double value = 0.0;
    std::for_each(p_literal_scalar_expression->begin(),
                  p_literal_scalar_expression->end(), [&value](auto& rValue) {
                      rValue = value;
                      value += 1.0;
                  });
    std::for_each(p_literal_array_3_expression->begin(),
                  p_literal_array_3_expression->end(), [&value](auto& rValue) {
                      rValue = value;
                      value += 1.0;
                  });

    auto p_binary_expression = BinaryExpression<BinaryOperations::Multiplication>::Create(p_literal_array_3_expression, p_literal_scalar_expression);

    std::size_t local_index = 0;
    std::for_each(p_binary_expression->cbegin(), p_binary_expression->cend(), [&p_binary_expression, &local_index](const auto Value) {
        std::size_t entity_index = local_index / 3;
        std::size_t component_index = local_index % 3;
        std::size_t entity_data_begin_index = entity_index * 3;
        ++local_index;
        KRATOS_CHECK_NEAR(
            Value, p_binary_expression->Evaluate(entity_index, entity_data_begin_index, component_index),
            1e-9);
    });
}

KRATOS_TEST_CASE_IN_SUITE(ExpressionIteratorFlatChar, KratosCoreFastSuite) {
    auto p_literal_scalar_expression = LiteralFlatExpression<char>::Create(10, {});
    auto p_literal_array_3_expression = LiteralFlatExpression<char>::Create(10, {3});

    char value = 0.0;
    std::for_each(p_literal_scalar_expression->begin(), p_literal_scalar_expression->end(), [&value](auto& rValue) {
        rValue = value;
        value += 1;
    });

    auto scalar_begin = p_literal_scalar_expression->cbegin();
    for (std::size_t i = 0; i < 10; ++i) {
        KRATOS_CHECK_EQUAL(*(scalar_begin + i), static_cast<char>(i));
    }

    std::for_each(p_literal_array_3_expression->begin(), p_literal_array_3_expression->end(), [&value](auto& rValue) {
        rValue = value;
        value += 1;
    });

    auto vector_begin = p_literal_array_3_expression->cbegin();
    for (std::size_t i = 0; i < 10; ++i) {
        KRATOS_CHECK_EQUAL(*(vector_begin + i), static_cast<char>(i + 10));
    }
}

KRATOS_TEST_CASE_IN_SUITE(ExpressionIteratorFlatInt, KratosCoreFastSuite) {
    auto p_literal_scalar_expression = LiteralFlatExpression<int>::Create(10, {});
    auto p_literal_array_3_expression = LiteralFlatExpression<int>::Create(10, {3});

    int value = 0.0;
    std::for_each(p_literal_scalar_expression->begin(), p_literal_scalar_expression->end(), [&value](auto& rValue) {
        rValue = value;
        value += 1;
    });

    auto scalar_begin = p_literal_scalar_expression->cbegin();
    for (std::size_t i = 0; i < 10; ++i) {
        KRATOS_CHECK_EQUAL(*(scalar_begin + i), static_cast<int>(i));
    }

    std::for_each(p_literal_array_3_expression->begin(), p_literal_array_3_expression->end(), [&value](auto& rValue) {
        rValue = value;
        value += 1;
    });

    auto vector_begin = p_literal_array_3_expression->cbegin();
    for (std::size_t i = 0; i < 10; ++i) {
        KRATOS_CHECK_EQUAL(*(vector_begin + i), static_cast<int>(i + 10));
    }
}

KRATOS_TEST_CASE_IN_SUITE(ExpressionIteratorFlatDouble, KratosCoreFastSuite) {
    auto p_literal_scalar_expression = LiteralFlatExpression<double>::Create(10, {});
    auto p_literal_array_3_expression = LiteralFlatExpression<double>::Create(10, {3});

    double value = 0.0;
    std::for_each(p_literal_scalar_expression->begin(), p_literal_scalar_expression->end(), [&value](auto& rValue) {
        rValue = value;
        value += 1;
    });

    auto scalar_begin = p_literal_scalar_expression->cbegin();
    for (std::size_t i = 0; i < 10; ++i) {
        KRATOS_CHECK_EQUAL(*(scalar_begin + i), static_cast<double>(i));
    }

    std::for_each(p_literal_array_3_expression->begin(), p_literal_array_3_expression->end(), [&value](auto& rValue) {
        rValue = value;
        value += 1;
    });

    auto vector_begin = p_literal_array_3_expression->cbegin();
    for (std::size_t i = 0; i < 10; ++i) {
        KRATOS_CHECK_EQUAL(*(vector_begin + i), static_cast<double>(i + 10));
    }
}

} // namespace Testing.
} // namespace Kratos.
