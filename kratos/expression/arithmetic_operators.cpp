//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

// Project includes
#include "expression/arithmetic_operators.h"
#include "expression/literal_expression.h"
#include "expression/binary_expression.h"


namespace Kratos {


#define KRATOS_DEFINE_BINARY_EXPRESSION_OPERATOR(OPERATOR_NAME, OPERATOR_CLASS)              \
    Expression::Pointer OPERATOR_NAME(const Expression::ConstPointer& rpLeft, const double Right) \
    {                                                                                        \
        return BinaryExpression<OPERATOR_CLASS>::Create(                                     \
            rpLeft, LiteralExpression<double>::Create(Right, rpLeft->NumberOfEntities()));   \
    }                                                                                        \
                                                                                             \
    Expression::Pointer OPERATOR_NAME(const double Left, const Expression::ConstPointer& rpRight) \
    {                                                                                        \
        return BinaryExpression<OPERATOR_CLASS>::Create(                                     \
            LiteralExpression<double>::Create(Left, rpRight->NumberOfEntities()), rpRight);  \
    }                                                                                        \
                                                                                             \
    Expression::Pointer OPERATOR_NAME(const Expression::ConstPointer& rpLeft,                     \
                                      const Expression::ConstPointer& rpRight)                    \
    {                                                                                        \
        KRATOS_ERROR_IF_NOT(                                                                 \
            rpLeft->NumberOfEntities() * rpLeft->GetItemComponentCount() ==                  \
            rpRight->NumberOfEntities() * rpRight->GetItemComponentCount())                  \
            << "Operand size mismatch in binary operator: " << #OPERATOR_NAME << "!\n"       \
            << "Left operand: " << *rpLeft << '\n'                                           \
            << "Right operand: " << *rpRight;                                                \
        return BinaryExpression<OPERATOR_CLASS>::Create(rpLeft, rpRight);                    \
    }

KRATOS_DEFINE_BINARY_EXPRESSION_OPERATOR(operator+, BinaryOperations::Addition)

KRATOS_DEFINE_BINARY_EXPRESSION_OPERATOR(operator-, BinaryOperations::Substraction)

KRATOS_DEFINE_BINARY_EXPRESSION_OPERATOR(operator*, BinaryOperations::Multiplication)

KRATOS_DEFINE_BINARY_EXPRESSION_OPERATOR(operator/, BinaryOperations::Division)

#undef KRATOS_DEFINE_BINARY_EXPRESSION_OPERATOR

} // namespace Kratos
